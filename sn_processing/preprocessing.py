from pathlib import Path
from specutils import Spectrum  # instead of Spectrum1D
import numpy as np
import astropy.units as u
from specutils.manipulation import FluxConservingResampler
from specutils.analysis import template_redshift  # <-- Use built-in function
from scipy.ndimage import convolve1d
import h5py
from typing import Tuple
import warnings
from astropy.utils.exceptions import AstropyWarning
import re
import matplotlib.pyplot as plt

warnings.simplefilter('ignore', AstropyWarning)

def read_spectrum(filepath: Path) -> Spectrum:
    """Reads a spectrum file into a Spectrum object."""
    if filepath.suffix == ".dat":
        # Read file, replace D/d with E for exponent notation
        with open(filepath, 'r') as f:
            lines = [line.replace('D', 'E').replace('d', 'E') for line in f]
        from io import StringIO
        data = np.loadtxt(StringIO(''.join(lines)))
        wavelength = data[:, 0] * u.AA
        flux = data[:, 1] * u.Unit("adu")  # Replace "adu" with correct units if known
        return Spectrum(spectral_axis=wavelength, flux=flux)
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
                flux = data * u.Unit("adu")  # Replace with correct flux unit if known

                return Spectrum(spectral_axis=wavelength, flux=flux)
    return None

def plot_spectra(spectrum_list):
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Resample each spectrum to match a common grid
    spectral_axis_values = [spec.spectral_axis.value for spec in spectrum_list]
    spectral_axis_min = min(min(spec) for spec in spectral_axis_values)
    spectral_axis_max = max(max(spec) for spec in spectral_axis_values)
    
    spectral_axis_common = np.linspace(spectral_axis_min, spectral_axis_max)
    
    flux_values_common = []
    for spec in spectrum_list:
        flux_values = np.interp(spectral_axis_common, spec.spectral_axis.value, spec.flux.value)
        flux_values_common.append(flux_values / np.max(spec.flux.value))
    
    # Plot each spectrum
    for i, flux_values in enumerate(flux_values_common):
        ax.plot(spectral_axis_common, flux_values + (i * 1), label=f"Spectrum {i+1}")
    
    ax.set_title("Spectra Plot")
    ax.set_xlabel("Wavelength (Å)")
    ax.set_ylabel("Normalized Flux")
    ax.legend()
    plt.show()

def find_best_redshift(observed_spec: Spectrum, template_spec: Spectrum, z_search_range: np.ndarray) -> float:
    """
    Finds the best redshift by cross-correlating an observed spectrum with a template,
    using specutils.analysis.template_redshift.
    """
    result = template_redshift(
        observed_spec,
        template_spec,
        z_search_range,
        resample_method='flux_conserving',
        extrapolation_treatment='truncate'
    )
    # result[1] is the best-fit redshift
    return result[1]

def bin_spectrum_constant_wavelength(spec: Spectrum, bin_width: u.Quantity) -> Spectrum:
    """Resamples a spectrum to a new grid with constant wavelength steps."""
    min_wl = spec.spectral_axis.min()
    max_wl = spec.spectral_axis.max()
    n_bins = int(((max_wl - min_wl) / bin_width).to_value(u.dimensionless_unscaled))
    new_grid = np.linspace(min_wl.value, max_wl.value, n_bins) * min_wl.unit
    resampler = FluxConservingResampler(extrapolation_treatment='zero_fill')
    return resampler(spec, new_grid)

def pad_spectrum_for_wavelet(spec: Spectrum, final_size: int = 4096) -> Tuple[np.ndarray, np.ndarray]:
    """Pads the spectrum to a power of two, matching IDL logic."""
    current_size = len(spec.flux)
    if current_size >= final_size:
        return spec.spectral_axis.value, spec.flux.value

    padded_flux = np.zeros(final_size)
    padded_flux[:current_size] = spec.flux.value

    # Extrapolate wavelength axis
    delta_lam = (spec.spectral_axis[-1] - spec.spectral_axis[0]) / (current_size - 1)
    start_lam = spec.spectral_axis[0].value
    padded_wl = start_lam + np.arange(final_size) * delta_lam.value

    # Fill extra points with a decaying exponential to minimize edge effects
    # This is a simple approximation of the IDL logic.
    fill_value = np.mean(spec.flux.value[-10:]) # Mean of last 10 points
    decay_factor = 0.99
    for i in range(current_size, final_size):
        fill_value *= decay_factor
        padded_flux[i] = fill_value

    return padded_wl, padded_flux


def atrous_transform(signal: np.ndarray, num_scales: int) -> np.ndarray:
    """
    Performs a 1D à trous wavelet transform on a signal.

    Args:
        signal: The 1D input signal (e.g., the flux array).
        num_scales: The number of wavelet scales to compute.

    Returns:
        A 2D numpy array where rows are the wavelet coefficients for each
        scale, plus the final smoothed residual. Shape is (num_scales + 1, len(signal)).
    """
    # The B3 spline scaling function is used as the convolution kernel.
    kernel = np.array([1/16, 1/4, 3/8, 1/4, 1/16])
    wavelet_coeffs = []
    c_prev = signal.copy()

    for j in range(num_scales):
        # At each scale 'j', the kernel is upsampled by inserting 2**j - 1 zeros
        # between its elements. This is equivalent to dilating the kernel.
        stride = 2**j
        upsampled_kernel = np.zeros((len(kernel) - 1) * stride + 1)
        upsampled_kernel[::stride] = kernel

        # Convolve with the upsampled kernel. 'mirror' mode handles boundaries.
        c_next = convolve1d(c_prev, upsampled_kernel, mode='mirror')

        # The wavelet coefficients for this scale are the difference plane.
        w_j = c_prev - c_next
        wavelet_coeffs.append(w_j)

        c_prev = c_next

    # The final smoothed (residual) plane is also part of the output.
    wavelet_coeffs.append(c_prev)
    return np.array(wavelet_coeffs)

def calculate_chi_index(
    wavelet_sum: np.ndarray,
    spectral_axis: np.ndarray,
    wl_min: float,
    wl_max: float
) -> float:
    """
    Calculates the 'chi' spectral index for a feature in a given wavelength range.
    The definition is based on the notes for sn1a_template_20110517.pro.
    """
    # Find indices corresponding to the wavelength region
    feature_indices = np.where((spectral_axis >= wl_min) & (spectral_axis <= wl_max))[0]

    if len(feature_indices) == 0:
        return 0.0

    # Extract the wavelet sum in that region
    feature_wavelet_sum = wavelet_sum[feature_indices]

    # Calculate standard deviation in the region
    std_dev = np.std(feature_wavelet_sum)

    if std_dev == 0:
        return 0.0

    # The index is the sum of the absolute values, normalized by std dev.
    # This is a plausible interpretation of the IDL code's intent.
    # The exact definition might need verification from the original paper.
    chi_index = np.sum(np.abs(feature_wavelet_sum)) / std_dev
    return chi_index


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
