from pathlib import Path
from specutils import Spectrum  # instead of Spectrum1D
import numpy as np
import astropy.units as u
from specutils.manipulation import FluxConservingResampler
from scipy import signal
from scipy.interpolate import interp1d
import h5py
from typing import Tuple, List, Dict
import warnings
from astropy.utils.exceptions import AstropyWarning
import re
import matplotlib.pyplot as plt
from sn_processing import io
import extinction

warnings.simplefilter('ignore', AstropyWarning)

def read_spectrum(filepath: Path, epoch: float = None) -> Spectrum:
    """
    Reads a spectrum file into a specutils.Spectrum object, including uncertainty if available.

    If reading fails and an epoch is provided, it attempts to find and load a template spectrum for that epoch.
    """
    try:
        if filepath.suffix == ".dat":
            # Read file, replace D/d with E for exponent notation
            with open(filepath, 'r') as f:
                lines = [line.replace('D', 'E').replace('d', 'E') for line in f]
            from io import StringIO
            data = np.loadtxt(StringIO(''.join(lines)))
            wavelength = data[:, 0] * u.AA
            flux = data[:, 1] * u.Unit("adu")  # Placeholder unit
            uncertainty = None
            if data.shape[1] > 2:
                from astropy.nddata import StdDevUncertainty
                print(f"  Found 3 columns in .dat file, assuming 3rd is uncertainty.")
                uncertainty = StdDevUncertainty(data[:, 2] * flux.unit)
            return Spectrum(spectral_axis=wavelength, flux=flux, uncertainty=uncertainty)
        elif filepath.suffix in [".fits", ".fit", ".flm"]:
            from astropy.io import fits
            try:
                 # Try reading the spectrum normally
                spec = Spectrum.read(filepath)
                # Try to access the unit; if it fails, patch it
                try:
                    _ = spec.spectral_axis.unit.to(u.AA)
                    return spec
                except Exception:
                    pass  # Will patch below if needed
            except Exception:
                # If Spectrum.read fails, manually extract data and patch units
                with fits.open(filepath) as hdul:
                    hdr = hdul[0].header
                    from astropy.nddata import StdDevUncertainty
                    data = hdul[0].data
                    # Print raw header for debugging
                    print(f"Raw CUNIT1: '{hdr.get('CUNIT1', '')}', Raw BUNIT1: '{hdr.get('BUNIT1', '')}'")
                    def clean_unit(unit_str):
                        # Remove all non-letter characters and lowercase
                        return re.sub(r'[^a-z]', '', unit_str.strip().lower())
                    cunit1 = clean_unit(hdr.get('CUNIT1', ''))
                    bunit1 = clean_unit(hdr.get('BUNIT1', ''))
                    print(f"Cleaned CUNIT1: {cunit1}, Cleaned BUNIT1: {bunit1}")
                    # Guess the unit
                    if 'micron' in cunit1 or 'micron' in bunit1:
                        unit = u.micron
                    elif 'angstrom' in cunit1 or 'angstroem' in cunit1 or 'angstrom' in bunit1 or 'angstroem' in bunit1:
                        unit = u.AA
                    else:
                        unit = u.AA  # Default/fallback

                    # Try to get the wavelength array
                    # This assumes a linear wavelength solution; adjust if your FITS files are different!
                    crval1 = hdr.get('CRVAL1', 0)  # starting wavelength
                    cdelt1 = hdr.get('CDELT1', 1)  # wavelength step
                    naxis1 = hdr.get('NAXIS1', len(data))
                    wavelength = crval1 + cdelt1 * np.arange(naxis1)
                    wavelength = wavelength * unit
                    flux = data * u.Unit("adu")  # Placeholder, should ideally come from BUNIT
                    uncertainty = None
                    # Check for uncertainty in a second HDU, a common convention
                    if len(hdul) > 1 and hdul[1].data is not None:
                        print(f"  Found uncertainty data in HDU 1.")
                        # Assuming the uncertainty has the same unit as the flux
                        uncertainty = StdDevUncertainty(hdul[1].data * flux.unit)

                    return Spectrum(spectral_axis=wavelength, flux=flux, uncertainty=uncertainty)
    except Exception as e:
        print(f"Warning: Could not read spectrum file {filepath.name}: {e}")
        if epoch is not None:
            print(f"Attempting to find a template for epoch {epoch:+.2f}...")
            matches = io.find_template_match(epoch)
            if matches:
                template_path = matches[0]
                print(f"  Found template: {template_path.name}. Loading it instead.")
                # Recursive call. Pass epoch=None to prevent an infinite loop if the template also fails to read.
                return read_spectrum(template_path, epoch=None)
            else:
                print(f"  No matching template found for epoch {epoch:+.2f}.")
            
    return None

def plot_spectra(
    spectrum_list: List[Spectrum],
    epoch_list: List[float],
    output_path: Path,
    title: str
):
    """
    Plots a list of spectra, normalized and offset, and saves the plot to a file.
    """
    if not spectrum_list:
        print("  - No spectra to plot. Skipping plot generation.")
        return

    fig, ax = plt.subplots(figsize=(12, 8))

    # Find a common wavelength range for all spectra
    spectral_axis_values = [spec.spectral_axis.value for spec in spectrum_list]
    common_min = min(s.min() for s in spectral_axis_values)
    common_max = max(s.max() for s in spectral_axis_values)
    common_grid = np.linspace(common_min, common_max, 2000)

    # Sort spectra by epoch for a clean, ordered plot
    sorted_data = sorted(zip(epoch_list, spectrum_list), key=lambda pair: pair[0])

    # Plot each spectrum
    for i, (epoch, spec) in enumerate(sorted_data):
        interpolated_flux = np.interp(common_grid, spec.spectral_axis.value, spec.flux.value)
        non_zero_flux = interpolated_flux[interpolated_flux > 0]
        norm_factor = np.median(non_zero_flux) if len(non_zero_flux) > 0 else 1.0
        if norm_factor <= 0: norm_factor = 1.0
        normalized_flux = interpolated_flux / norm_factor
        offset = i * 0.8
        ax.plot(common_grid, normalized_flux + offset, label=f"Epoch {epoch:+.1f}")

    ax.set_title(title)
    ax.set_xlabel("Rest Wavelength (Å)")
    ax.set_ylabel("Normalized Flux + Offset")
    ax.legend(title="Epoch (days)", fontsize='small')

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  - Saved summary plot to {output_path}")


def measure_redshift_from_template(
    observed_spectrum: Spectrum,
    template_spectrum: Spectrum,
    resampling_points: int = 2000
) -> Tuple[float, float]:
    """
    Measures the redshift of an observed spectrum by cross-correlating with a template.

    The method works by:
    1. Finding the overlapping wavelength range between the two spectra.
    2. Resampling both spectra onto a common, linearly-spaced log-wavelength grid.
       A multiplicative shift in wavelength (due to redshift) becomes an additive
       shift in log-wavelength space, which is ideal for cross-correlation.
    3. Performing a cross-correlation on the resampled flux values.
    4. Finding the peak of the cross-correlation function.
    5. Converting the lag at the peak back into a redshift value.

    Args:
        observed_spectrum (Spectrum): The observed supernova spectrum (`specutils.Spectrum`).
        template_spectrum (Spectrum): The template spectrum at z=0 (`specutils.Spectrum`).
        resampling_points (int): The number of points for the log-wavelength grid.

    Returns:
        A tuple containing:
        - measured_redshift (float): The calculated redshift.
        - correlation_peak (float): The peak value of the cross-correlation function,
                                    indicating the quality of the fit.
    """
    # 1. Determine the overlapping wavelength range.
    common_wave_min = max(
        observed_spectrum.spectral_axis.min(), template_spectrum.spectral_axis.min()
    )
    common_wave_max = min(
        observed_spectrum.spectral_axis.max(), template_spectrum.spectral_axis.max()
    )

    if common_wave_min >= common_wave_max:
        raise ValueError("Observed spectrum and template have no overlapping wavelength range.")

    # 2. Create a common log-wavelength grid and interpolate both spectra onto it.
    log_wave_grid = np.linspace(
        np.log10(common_wave_min.value),
        np.log10(common_wave_max.value),
        resampling_points
    )

    obs_interp = interp1d(
        np.log10(observed_spectrum.spectral_axis.value),
        observed_spectrum.flux.value,
        bounds_error=False,
        fill_value=0.0
    )
    template_interp = interp1d(
        np.log10(template_spectrum.spectral_axis.value),
        template_spectrum.flux.value,
        bounds_error=False,
        fill_value=0.0
    )

    resampled_obs_flux = obs_interp(log_wave_grid)
    resampled_template_flux = template_interp(log_wave_grid)

    # Normalize fluxes for robust correlation
    resampled_obs_flux = (resampled_obs_flux - np.mean(resampled_obs_flux)) / np.std(resampled_obs_flux)
    resampled_template_flux = (resampled_template_flux - np.mean(resampled_template_flux)) / np.std(resampled_template_flux)

    # 3. Perform cross-correlation
    correlation = signal.correlate(resampled_obs_flux, resampled_template_flux, mode='same')

    # 4. Find the peak of the correlation
    peak_index = np.argmax(correlation)
    correlation_peak = correlation[peak_index]

    # 5. Convert the lag (shift in array index) to a redshift
    center_index = len(correlation) // 2
    lag = peak_index - center_index
    log_wave_step = log_wave_grid[1] - log_wave_grid[0]
    log_z_shift = lag * log_wave_step  # This is log10(1+z)

    measured_redshift = 10**log_z_shift - 1

    return measured_redshift, correlation_peak


def find_best_redshift(
    observed_spectrum: Spectrum,
    template_spectra: List[Spectrum],
    **kwargs
) -> Tuple[float, float, Spectrum]:
    """
    Finds the best-fit redshift by comparing an observed spectrum against multiple templates.

    This function iterates through a list of template spectra (e.g., representing
    different supernova epochs), measures the redshift against each one using
    cross-correlation, and selects the template that yields the highest
    correlation score.

    Args:
        observed_spectrum (Spectrum): The observed supernova spectrum (`specutils.Spectrum`).
        template_spectra (List[Spectrum]): A list of template spectra to match against.
        **kwargs: Keyword arguments to pass to `measure_redshift_from_template`.

    Returns:
        A tuple containing:
        - best_redshift (float): The redshift from the best-fitting template.
        - best_correlation_score (float): The peak correlation score for the best fit.
        - best_template (Spectrum): The template spectrum that gave the best fit.
    """
    if not template_spectra:
        raise ValueError("template_spectra list cannot be empty.")

    results = []
    for template in template_spectra:
        try:
            redshift, score = measure_redshift_from_template(observed_spectrum, template, **kwargs)
            results.append((score, redshift, template))
        except ValueError:
            continue  # Skip templates with no wavelength overlap

    if not results:
        raise RuntimeError("Could not find a suitable redshift for any of the provided templates.")

    # Find the result with the maximum correlation score
    best_score, best_redshift, best_template = max(results, key=lambda item: item[0])

    return best_redshift, best_score, best_template


def generate_hsiao_templates(dm15: float, hsiao_data: Dict[str, np.ndarray]) -> List[Spectrum]:
    """
    Generates a list of 2D spectra from the 3D Hsiao template for a specific dm15.

    This function interpolates the Hsiao flux cube at the given dm15 value to
    produce a set of template spectra, one for each phase in the template grid.

    Args:
        dm15: The target dm15 value to interpolate at.
        hsiao_data: The loaded Hsiao data cube from `io.load_hsiao_data_cube`.

    Returns:
        A list of `specutils.Spectrum` objects, one for each phase.
    """
    if not hsiao_data:
        return []

    wavelengths = hsiao_data['wavelengths']
    phases = hsiao_data['phases']
    dm15s = hsiao_data['dm15s']
    flux_cube = hsiao_data['flux_cube']  # Shape: (n_phases, n_dm15s, n_wavelengths)

    generated_templates = []
    # For each phase, we have a 2D slice of (dm15 vs. wavelength)
    for i, phase in enumerate(phases):
        flux_slice_2d = flux_cube[i, :, :]  # Shape: (n_dm15s, n_wavelengths)

        # Create an interpolation function for this phase slice.
        # It will interpolate along the dm15 axis (axis=0).
        interp_func = interp1d(
            dm15s, flux_slice_2d, axis=0, bounds_error=False, fill_value='extrapolate'
        )

        # Get the interpolated flux array for the target dm15
        interpolated_flux = interp_func(dm15)

        # Create a Spectrum object (assuming arbitrary flux units for the template)
        template_spec = Spectrum(spectral_axis=wavelengths*u.AA, flux=interpolated_flux*u.Unit('adu'))
        generated_templates.append(template_spec)

    return generated_templates

def correct_for_redshift(spectrum: Spectrum, z: float) -> Spectrum:
    """
    Corrects a spectrum to its rest frame (de-redshifts).

    Args:
        spectrum: The observed spectrum.
        z: The redshift of the source.

    Returns:
        A new Spectrum object in the rest frame.
    """
    if z < 0:
        raise ValueError("Redshift z must be non-negative.")
    if z == 0:
        return spectrum

    # De-redshift the wavelength axis: lambda_rest = lambda_obs / (1 + z)
    rest_wavelength = spectrum.spectral_axis / (1 + z)

    # Correct flux for cosmological expansion: F_rest = F_obs * (1 + z)
    rest_flux = spectrum.flux * (1 + z)

    # Propagate uncertainty if it exists
    rest_uncertainty = None
    if spectrum.uncertainty is not None:
        from astropy.nddata import StdDevUncertainty
        rest_uncertainty_array = spectrum.uncertainty.array * (1 + z)
        rest_uncertainty = StdDevUncertainty(rest_uncertainty_array * rest_flux.unit)

    return Spectrum(
        spectral_axis=rest_wavelength,
        flux=rest_flux,
        uncertainty=rest_uncertainty
    )


def correct_for_reddening(spectrum: Spectrum, ebv: float, r_v: float = 3.1) -> Spectrum:
    """
    Corrects a spectrum for interstellar reddening using the CCM89 law.

    Args:
        spectrum: The observed, reddened spectrum.
        ebv: E(B-V) color excess.
        r_v: Ratio of total to selective extinction.

    Returns:
        A new, de-reddened Spectrum object.
    """
    if ebv < 0:
        raise ValueError("E(B-V) must be non-negative.")
    if ebv == 0:
        return spectrum

    # Wavelength must be in Angstroms for the extinction package
    wave_angstrom = spectrum.spectral_axis.to_value(u.AA)

    # Calculate the extinction in magnitudes for each wavelength
    a_lambda = extinction.ccm89(wave_angstrom, a_v=ebv * r_v, r_v=r_v)

    # Apply the de-reddening correction to the flux
    dereddened_flux_values = extinction.remove(a_lambda, spectrum.flux.value)
    dereddened_flux = dereddened_flux_values * spectrum.flux.unit

    # Propagate uncertainty. De-reddening is a multiplicative factor.
    dereddened_uncertainty = None
    if spectrum.uncertainty is not None:
        from astropy.nddata import StdDevUncertainty
        # The correction factor is dereddened_flux / original_flux
        correction_factor = dereddened_flux_values / spectrum.flux.value
        # Replace NaNs or Infs from division by zero with 1 (no change)
        correction_factor[~np.isfinite(correction_factor)] = 1.0
        dereddened_uncertainty_array = spectrum.uncertainty.array * correction_factor
        dereddened_uncertainty = StdDevUncertainty(dereddened_uncertainty_array * dereddened_flux.unit)

    return Spectrum(
        spectral_axis=spectrum.spectral_axis,
        flux=dereddened_flux,
        uncertainty=dereddened_uncertainty
    )


def bin_spectrum_constant_wavelength(spec: Spectrum, bin_width: u.Quantity) -> Spectrum:
    """Resamples a spectrum to a new grid with constant wavelength steps."""
    min_wl = spec.spectral_axis.min()
    max_wl = spec.spectral_axis.max()
    n_bins = int(((max_wl - min_wl) / bin_width).to_value(u.dimensionless_unscaled))
    new_grid = np.linspace(min_wl.value, max_wl.value, n_bins) * min_wl.unit
    resampler = FluxConservingResampler(extrapolation_treatment='zero_fill')
    return resampler(spec, new_grid)

def pad_spectrum_for_wavelet(spec: Spectrum, final_size: int = 4096) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Pads the spectrum's wavelength, flux, and uncertainty to a power of two.
    """
    current_size = len(spec.flux)
    if current_size >= final_size:
        uncertainty_array = spec.uncertainty.array if spec.uncertainty is not None else np.zeros_like(spec.flux.value)
        return spec.spectral_axis.value, spec.flux.value, uncertainty_array

    padded_flux = np.zeros(final_size)
    padded_uncertainty = np.zeros(final_size)
    padded_flux[:current_size] = spec.flux.value

    # Extrapolate wavelength axis
    delta_lam = (spec.spectral_axis[-1] - spec.spectral_axis[0]) / (current_size - 1)
    start_lam = spec.spectral_axis[0].value
    padded_wl = start_lam + np.arange(final_size) * delta_lam.value

    # Handle uncertainty padding
    if spec.uncertainty is not None:
        padded_uncertainty[:current_size] = spec.uncertainty.array
        # Fill the rest with a high value (max uncertainty) to indicate low confidence
        fill_value = np.max(spec.uncertainty.array)
        padded_uncertainty[current_size:] = fill_value
    else:
        # If no uncertainty, the padded array remains zeros
        pass

    # Fill extra points with a decaying exponential to minimize edge effects
    # This is a simple approximation of the IDL logic.
    flux_fill_value = np.mean(spec.flux.value[-10:]) # Mean of last 10 points
    decay_factor = 0.99
    for i in range(current_size, final_size):
        flux_fill_value *= decay_factor
        padded_flux[i] = flux_fill_value

    return padded_wl, padded_flux, padded_uncertainty


def add_poisson_noise(spectrum: Spectrum, target_snr: float) -> Spectrum:
    """
    Adds Poisson-like noise to a spectrum to achieve a target signal-to-noise ratio.

    This function simulates the noise profile of a photon-counting instrument.
    It scales the flux to an equivalent number of "counts" where the desired SNR
    is met at a reference flux level (here, the median flux), adds Poisson noise,
    and then scales the flux and new uncertainty back to the original units.

    Args:
        spectrum (Spectrum): The input, relatively noise-free spectrum.
        target_snr (float): The desired signal-to-noise ratio at the median flux level.

    Returns:
        A new Spectrum object with the added noise and a corresponding uncertainty array.
    """
    if target_snr <= 0:
        raise ValueError("Target SNR must be positive.")

    # Use the median flux as the reference signal level.
    # Ensure we don't take the median of zero-padded regions.
    non_zero_flux = spectrum.flux[spectrum.flux.value > 0]
    if len(non_zero_flux) == 0:
        # If all flux is zero, return the original spectrum
        return spectrum

    ref_flux_val = np.median(non_zero_flux.value)
    if ref_flux_val <= 0:
        # Handle cases where median is still zero or negative by falling back to mean
        ref_flux_val = np.mean(non_zero_flux.value)
        if ref_flux_val <= 0:
             raise ValueError("Cannot add noise to a spectrum with no positive flux.")

    # For Poisson noise, SNR = sqrt(signal_in_counts).
    # So, signal_in_counts = SNR^2.
    # We find a scaling factor 'k' to convert flux to counts such that
    # k * ref_flux_val = target_snr**2.
    k = target_snr**2 / ref_flux_val

    # Convert the entire flux array to "effective counts".
    # Ensure counts are non-negative for the Poisson distribution.
    flux_in_counts = np.maximum(0, spectrum.flux.value * k)

    # Generate a new flux array by sampling from a Poisson distribution.
    noisy_flux_in_counts = np.random.poisson(flux_in_counts)

    # Convert the noisy counts back to the original flux units.
    noisy_flux_values = noisy_flux_in_counts / k
    noisy_flux = noisy_flux_values * spectrum.flux.unit

    # The uncertainty of a Poisson process is sqrt(counts).
    # Propagate this back to flux units.
    from astropy.nddata import StdDevUncertainty
    # Use max(1,...) to avoid sqrt(0) which can be problematic.
    uncertainty_values = np.sqrt(np.maximum(1, noisy_flux_in_counts)) / k
    uncertainty = StdDevUncertainty(uncertainty_values * spectrum.flux.unit)

    # Create and return the new spectrum object.
    return Spectrum(
        spectral_axis=spectrum.spectral_axis,
        flux=noisy_flux,
        uncertainty=uncertainty
    )


def generate_monte_carlo_spectra(
    base_spectrum: Spectrum, snr_levels: List[float], num_realizations: int = 1
) -> Dict[float, List[Spectrum]]:
    """
    Generates multiple realizations of a spectrum with different noise levels for Monte Carlo.
    """
    mc_spectra = {}
    for snr in snr_levels:
        mc_spectra[snr] = []
        for _ in range(num_realizations):
            noisy_spec = add_poisson_noise(base_spectrum, snr)
            mc_spectra[snr].append(noisy_spec)
    return mc_spectra


def save_results(output_file: str, sn_name: str, epoch: float, data_dict: dict):
    """
    Saves processed data for a single epoch to an HDF5 file.

    Args:
        output_file: Path to the HDF5 file.
        sn_name: Name of the supernova.
        epoch: The observation epoch.
        data_dict: A dictionary of numpy arrays to save (e.g., {'wavelets': w, 'chi': c}).
    """
    with h5py.File(output_file, 'a') as f:  # 'a' mode will create or append
        # Create a group for the supernova if it doesn't exist
        sn_group = f.require_group(sn_name)
        # Create a group for the epoch, overwriting if it exists
        epoch_str = f"epoch_{epoch:+.2f}"
        if epoch_str in sn_group:
            del sn_group[epoch_str]
        epoch_group = sn_group.create_group(epoch_str)

        # Save each item in the data_dict as a dataset
        for key, value in data_dict.items():
            epoch_group.create_dataset(key, data=value)

        print(f"Saved data for {sn_name} at epoch {epoch}")
