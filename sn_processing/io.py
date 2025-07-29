import re
from pathlib import Path
from typing import List, Dict
from datetime import datetime
import pandas as pd
from functools import lru_cache
from typing import Any
from astropy.io import fits

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
    # Define column names based on the dm15.dat file structure
    col_names = ['SN_name', 'SN_dm15', 'SN_dm15source', 'SN_datasource', 'SN_z']
    try:
        df = pd.read_csv(param_file,
                         delimiter=r'\s+',  # Use regex for one or more spaces
                         names=col_names,
                         header=None,
                         comment='#')
        # Clean and set the SN name as the index for easy lookup
        df['SN_name'] = df['SN_name'].str.strip().str.lower()
        df.set_index('SN_name', inplace=True)
        return df
    except FileNotFoundError:
        print(f"Warning: Parameter file not found at {param_file}")
        return pd.DataFrame()

def get_sn_redshift(sn_name: str) -> float:
    """Looks up the redshift for a given supernova from the parameter file."""
    params_df = load_sn_parameters()
    if params_df.empty:
        return None
    try:
        return float(params_df.loc[sn_name.lower(), 'SN_z'])
    except (KeyError, ValueError):
        print(f"Warning: Redshift not found for {sn_name}")
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
