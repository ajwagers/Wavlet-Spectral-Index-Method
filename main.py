from sn_processing import io, preprocessing, wavelets
import astropy.units as u
import numpy as np
from pathlib import Path

# --- Configuration ---
DATA_DIR = "./SNSPEC"
OUTPUT_FILE = "./results/spectral_features.h5"

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
NUM_WAVELET_SCALES = 5 # Example value
SCALES_TO_SUM = [2, 3, 4] # Scales to sum for feature extraction

# This would typically be loaded from a file (e.g., a CSV or JSON).
# See sn_processing/io.py for more details.
# Dates are in Julian Day format.
SUPERNOVA_MAX_LIGHT_DATES = {
    'sn2001el': 2452180.5, # Corresponds to roughly 2001-09-28
    # Add other supernovae and their max light dates (in JD) here
}

def main():
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

            # Add the current spectrum to the list
            spectrum_list.append(spec)

            # 2. TODO: Redshift correction (if needed)
            # spec = preprocessing.correct_redshift(spec, ...)

            # 3. Bin and Pad
            # Example: bin to 5 Angstroms
            binned_spec = preprocessing.bin_spectrum_constant_wavelength(spec, 5 * u.AA)
            if binned_spec.uncertainty is not None:
                print(f"    - Propagating uncertainties through binning and padding.")
            else:
                print(f"    - No uncertainty data found to propagate.")
            wl, flux, uncertainty = preprocessing.pad_spectrum_for_wavelet(binned_spec)

            # 4. Wavelet Transform
            wavelet_coeffs = wavelets.atrous_transform(flux, NUM_WAVELET_SCALES)

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
