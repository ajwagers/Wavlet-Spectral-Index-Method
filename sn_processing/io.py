import re
from pathlib import Path
from typing import List, Dict
from datetime import datetime
import pandas as pd
from functools import lru_cache
from typing import Any
from astropy.io import fits
import numpy as np

def find_template_match(epoch: float) -> List[Path]:
    """Finds template spectra in various libraries that match the given epoch."""
    template_libraries = {
        "HSIAO": Path("/SNSPEC/HSIAO/"),
        "MODEL": Path("/SNSPEC/MODEL/"),
        "NUGENT": Path("/SNSPEC/NUGENT/"),
    }

    all_matches = []
    for model_name, template_dir in template_libraries.items():
        if not template_dir.is_dir():
            print(f"Warning: Template directory for {model_name} not found: {template_dir}")
            continue

        for ext in ("dat", "flm", "fits", "fit"):
            for fpath in template_dir.glob(f"*.{ext}"):
                # Try to parse epoch from filename, e.g., ".p015" or ".m04"
                match = re.search(r'\.(m|p)(\d+)', fpath.stem)
                if not match:
                    continue

                sign = -1.0 if match.group(1) == 'm' else 1.0
                file_epoch = sign * int(match.group(2))

                if file_epoch == epoch:
                    all_matches.append(fpath)

    return all_matches

def find_spectra_for_sn(base_dir: str, sn_name: str) -> List[Path]:
    """Finds all spectrum files for a given supernova."""
    folder_name = sn_name[2:] if sn_name.lower().startswith("sn") else sn_name
    search_path = Path(base_dir) / folder_name
    if not search_path.exists():
        print(f"Warning: Directory not found for {sn_name} at {search_path}")
        return []
    # Lowercase for matching
    sn_name_lower = sn_name.lower()
    # Use a set to automatically handle duplicates
    found_files = set()
    for ext in ("dat", "flm", "fits", "fit"):
        found_files.update(search_path.glob(f"{sn_name_lower}.*.{ext}"))

    # The original code passed a compiled regex to glob, causing a TypeError.
    # The regex was equivalent to matching any .dat file, which we do here.
    # This is likely to find spectra with different naming conventions.
    found_files.update(search_path.glob('*.dat'))

    # Filter out hidden files (e.g., .DS_Store on macOS) that start with a dot.
    visible_files = [f for f in found_files if not f.name.startswith('.')]
    return sorted(visible_files)

def parse_epoch_from_filename(filepath: Path) -> float:
    """
    Extracts the epoch (days from max light) from a filename like sn2001el.m04.dat or sn2001el.p01.dat.
    Returns the epoch as a float (negative for 'm', positive for 'p').
    """
    import re
    match = re.search(r'\.(m|p)(\d+)', filepath.stem)
    if not match:
        raise ValueError(f"Could not parse epoch from filename: {filepath.name}")

    sign = -1 if match.group(1) == 'm' else 1
    epoch = sign * int(match.group(2))
    return float(epoch)

@lru_cache(maxsize=1)
def load_sn_parameters(param_file: Path = Path('./dm15.dat')) -> pd.DataFrame:
    """
    Loads supernova parameter data from a file into a DataFrame, caching the result.
    The file is expected to be a comma-delimited file without a header.
    """
    # Define a comprehensive list of column names for the dm15.dat file.
    # Pandas will handle rows that have fewer columns gracefully.
    col_names = [
        'SN_name', 'dm15', 'dm15_source', 'data_source', 'z', 'z_source',
        'Bmax', 'Bmax_err', 'Bmax_source', 'Vmax', 'Vmax_err', 'Vmax_source',
        'E(B-V)_MW', 'Host_morphology', 'E(B-V)_host', 'E(B-V)_host_err',
        'E(B-V)_host_source', 'NaI_D_EW', 'NaI_D_EW_err', 'NaI_D_EW_source',
        'SiII_vel', 'SiII_vel_grad', 'SiII_EW_1', 'SiII_EW_2', 'SiII_EW_3', 'SiII_source'
    ]
    try:
        df = pd.read_csv(param_file,
                         delimiter=',',      # The file is comma-delimited
                         names=col_names,
                         header=None,
                         comment='#')
        # Clean up the data
        df['SN_name'] = df['SN_name'].str.strip().str.lower()
        # Convert relevant columns to numeric, coercing errors to NaN
        df['z'] = pd.to_numeric(df['z'], errors='coerce')
        df['dm15'] = pd.to_numeric(df['dm15'], errors='coerce')
        df.set_index('SN_name', inplace=True)
        return df
    except FileNotFoundError:
        print(f"Warning: Parameter file not found at {param_file}")
        return pd.DataFrame()

@lru_cache(maxsize=1)
def load_hsiao_data_cube(template_path: Path) -> Dict[str, np.ndarray]:
    """
    Loads the 3D Hsiao template data from a NumPy .npz file.

    The .npz file should contain the following arrays:
    - 'wavelengths': 1D array of wavelength points.
    - 'phases': 1D array of phase values (days from max light).
    - 'dm15s': 1D array of dm15 values.
    - 'flux_cube': 3D array of flux with shape (n_phases, n_dm15s, n_wavelengths).

    Returns:
        A dictionary containing the template data, or an empty dictionary if not found.
    """
    if not template_path.exists():
        print(f"Warning: Hsiao template data cube not found at {template_path}")
        return {}
    try:
        data = np.load(template_path)
        # Check for required keys
        required_keys = ['wavelengths', 'phases', 'dm15s', 'flux_cube']
        if not all(key in data for key in required_keys):
            raise KeyError(f"Hsiao template file is missing one of the required keys: {required_keys}")
        return data
    except Exception as e:
        print(f"Error loading Hsiao template data cube: {e}")
        return {}

def get_sn_redshift(sn_name: str) -> float:
    """Looks up the redshift for a given supernova from the parameter file."""
    params_df = load_sn_parameters()
    if params_df.empty:
        return None
    try:
        redshift = params_df.loc[sn_name.lower(), 'z']
        return float(redshift) if pd.notna(redshift) else None
    except (KeyError, ValueError):
        print(f"Warning: Redshift not found for {sn_name}")
        return None

def get_sn_dm15(sn_name: str) -> float:
    """Looks up the dm15 value for a given supernova from the parameter file."""
    params_df = load_sn_parameters()
    if params_df.empty:
        return None
    try:
        dm15 = params_df.loc[sn_name.lower(), 'dm15']
        return float(dm15) if pd.notna(dm15) else None
    except (KeyError, ValueError):
        return None

def extract_fits_metadata(filepath: Path) -> Dict[str, Any]:
    """
    Extracts key metadata from a FITS file header.

    This function opens a FITS file, reads the primary header, and extracts
    a predefined set of common metadata keywords like object name and
    observation date.

    Args:
        filepath: The path to the FITS file.

    Returns:
        A dictionary containing the metadata found. Returns an empty dictionary
        if the file is not a FITS file or the header cannot be read.
    """
    if filepath.suffix.lower() not in ['.fits', '.fit', '.flm']:
        return {}

    try:
        with fits.open(filepath) as hdul:
            header = hdul[0].header
            metadata = {
                'object': header.get('OBJECT'),
                'date_obs': header.get('DATE-OBS'),
                'exptime': header.get('EXPTIME'),
            }
            # Filter out keys where the value is None
            return {k: v for k, v in metadata.items() if v is not None}
    except Exception as e:
        print(f"Warning: Could not read FITS header from {filepath.name}: {e}")
        return {}
