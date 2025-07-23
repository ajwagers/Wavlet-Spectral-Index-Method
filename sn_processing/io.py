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
    # Match .dat, .flm, .fits, .fit (dot after SN name)
    files = []
    for ext in ("dat", "flm", "fits", "fit"):
        files += list(search_path.glob(f"{sn_name_lower}.*.{ext}"))
    return sorted(files)

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
