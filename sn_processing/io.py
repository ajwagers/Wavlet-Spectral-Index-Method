import re
from pathlib import Path
from typing import List, Dict
from datetime import datetime

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
