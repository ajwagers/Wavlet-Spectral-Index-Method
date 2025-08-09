from sn_processing import io, preprocessing, wavelets
import astropy.units as u
import numpy as np
from pathlib import Path

# --- Configuration ---
DATA_DIR = "./SNSPEC"
OUTPUT_FILE = "./results/spectral_features.h5"
WAVELET_PLOT_DIR = "./results/wavelet_plots"

# Automatically detect all supernovae folders in SNSPEC
SUPERNOVAE_TO_PROCESS = [
#    f"sn{folder.name}" for folder in Path(DATA_DIR).iterdir()
#    if folder.is_dir() and not folder.name.startswith('.')
"sn1981B","sn1983G","sn1984A","sn1986G","sn1989B","sn1990N","sn1990O","sn1991bg","sn2001el"
]

FEATURE_DEFINITIONS = {
    'SiII_6355': (6000, 6250), # Example wavelength ranges in Angstrom
    'SiII_5485': (5300, 5500),
    # Add other features from your notes...
}
NUM_WAVELET_SCALES = 10 # Match the IDL implementation
SCALES_TO_SUM = [2, 3, 4] # Scales to sum for feature extraction

# This would typically be loaded from a file (e.g., a CSV or JSON).
# See sn_processing/io.py for more details.
# Dates are in Julian Day format.
SUPERNOVA_MAX_LIGHT_DATES = {
    'sn2001el': 2452180.5, # Corresponds to roughly 2001-09-28
    # Add other supernovae and their max light dates (in JD) here
}

# This would also typically be loaded from a parameter file.
# For now, we define it here. These are example values.
SUPERNOVA_EBV_VALUES = {
    'sn1981b': 0.12,
    'sn1989b': 0.35,
    'sn1990n': 0.0,
    'sn1991bg': 0.1,
    'sn2001el': 0.068,
}

def main():
    # Create the main plot directory to store visualizations
    Path(WAVELET_PLOT_DIR).mkdir(parents=True, exist_ok=True)

    for sn_name in SUPERNOVAE_TO_PROCESS:
        print(f"--- Processing {sn_name} ---")
        spectrum_files = io.find_spectra_for_sn(DATA_DIR, sn_name)

        if not spectrum_files:
            print(f"No spectra found for {sn_name}. Skipping.")
            continue

        spectrum_list = []
        for spec_file in spectrum_files:
            try:
                # Pass the lookup dictionary to the parsing function
                epoch = io.parse_epoch_from_filename(spec_file)
            except ValueError as e:
                print(f"  Skipping {spec_file.name}: {e}")
                continue
            print(f"  Processing epoch: {epoch:+.2f}")

            # 1. Read spectrum and metadata
            metadata = io.extract_fits_metadata(spec_file)
            if metadata.get('object'):
                print(f"    - Object from FITS header: {metadata['object']}")
            spec = preprocessing.read_spectrum(spec_file, epoch=epoch)
            if spec is None:
                continue

            # 2. Apply Corrections (De-redshifting and De-reddening)
            # The correct order is to first correct for redshift to get to the
            # rest frame, then correct for reddening in the rest frame.
            redshift = io.get_sn_redshift(sn_name)
            ebv = SUPERNOVA_EBV_VALUES.get(sn_name.lower(), 0.0)

            corrected_spec = spec
            if redshift is not None and redshift > 0:
                print(f"    - Correcting for redshift z = {redshift:.4f}")
                corrected_spec = preprocessing.correct_for_redshift(corrected_spec, z=redshift)

            if ebv > 0:
                print(f"    - Correcting for reddening with E(B-V) = {ebv:.3f}")
                corrected_spec = preprocessing.correct_for_reddening(corrected_spec, ebv=ebv)

            spectrum_list.append(corrected_spec)

            # 3. Bin and Pad
            # Use the corrected spectrum for the rest of the analysis
            binned_spec = preprocessing.bin_spectrum_constant_wavelength(corrected_spec, 5 * u.AA)
            if binned_spec.uncertainty is not None:
                print(f"    - Propagating uncertainties through binning and padding.")
            else:
                print(f"    - No uncertainty data found to propagate.")
            wl, flux, uncertainty = preprocessing.pad_spectrum_for_wavelet(binned_spec)

            # 4. Wavelet Transform
            wavelet_coeffs = wavelets.atrous_transform(flux, NUM_WAVELET_SCALES)

            # 4.5 Visualize Decomposition and save to file
            plot_filename = f"{sn_name}_epoch_{epoch:+.2f}.png"
            plot_output_path = Path(WAVELET_PLOT_DIR) / plot_filename
            wavelets.plot_wavelet_decomposition(
                spectral_axis=wl,
                original_signal=flux,
                wavelet_coeffs=wavelet_coeffs,
                output_path=plot_output_path,
                title=f"Wavelet Decomposition for {sn_name} at Epoch {epoch:+.2f}"
            )

            # 5. Feature Measurement
            # The IDL code summed scales 2, 3, and 4 (0-indexed: 1, 2, 3)
            wavelet_sum = wavelets.sum_wavelet_scales(
                wavelet_coeffs, scales_to_sum=SCALES_TO_SUM
            )
            chi_indices = wavelets.extract_features(
                wavelet_sum, wl, FEATURE_DEFINITIONS
            )

            # 6. Save results
            results_to_save = {
                'wavelengths': wl,
                'binned_flux': flux,
                'binned_uncertainty': uncertainty,
                'wavelet_coeffs': wavelet_coeffs,
                'wavelet_sum_234': wavelet_sum,
                'chi_indices': np.array(list(chi_indices.values())),
                'chi_names': np.array(list(chi_indices.keys()), dtype='S')
            }
            preprocessing.save_results(OUTPUT_FILE, sn_name, epoch, results_to_save)
        
        # Now plot all spectra
        preprocessing.plot_spectra(spectrum_list)

if __name__ == "__main__":
    main()
