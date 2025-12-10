"""
DSSP Parser for ABC Domain Predictor

Parses DSSP output to extract secondary structure assignments and
hydrogen bond information for enhancing the contact graph.

DSSP format reference:
- Column 1-5: Residue number
- Column 6: Chain ID
- Column 7-10: Amino acid (3-letter code)
- Column 17: Secondary structure (H=helix, E=strand, T=turn, S=bend, G=310-helix, B=bridge, I=pi-helix, ' '=coil)
- Column 26-38: NH-->O hydrogen bond (residue offset, energy)
- Column 39-50: O-->HN hydrogen bond (residue offset, energy)
- Column 51-62: NH-->O hydrogen bond #2
- Column 63-74: O-->HN hydrogen bond #2
"""

import logging
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


@dataclass
class HydrogenBond:
    """A hydrogen bond between two residues."""
    donor_resnum: int      # Residue donating NH
    acceptor_resnum: int   # Residue accepting O
    energy: float          # H-bond energy (kcal/mol, more negative = stronger)

    @property
    def is_strong(self) -> bool:
        """Strong H-bond (energy < -0.5 kcal/mol)."""
        return self.energy < -0.5

    @property
    def sequence_separation(self) -> int:
        """Absolute sequence distance between bonded residues."""
        return abs(self.acceptor_resnum - self.donor_resnum)

    @property
    def is_long_range(self) -> bool:
        """Long-range H-bond (|i-j| > 5), typically beta sheet."""
        return self.sequence_separation > 5


@dataclass
class DSSPResidue:
    """DSSP information for a single residue."""
    resnum: int
    chain: str
    aa: str
    ss: str  # Secondary structure code
    acc: int  # Solvent accessibility
    hbonds: List[HydrogenBond]  # H-bonds involving this residue


@dataclass
class DSSPData:
    """Complete DSSP analysis results."""
    residues: Dict[int, DSSPResidue]  # resnum -> DSSPResidue
    hbonds: List[HydrogenBond]  # All hydrogen bonds

    @property
    def long_range_hbonds(self) -> List[HydrogenBond]:
        """Get only long-range H-bonds (beta sheet contacts)."""
        return [hb for hb in self.hbonds if hb.is_long_range]

    @property
    def strong_hbonds(self) -> List[HydrogenBond]:
        """Get only strong H-bonds."""
        return [hb for hb in self.hbonds if hb.is_strong]

    def get_ss(self, resnum: int) -> str:
        """Get secondary structure for a residue."""
        if resnum in self.residues:
            return self.residues[resnum].ss
        return ' '

    def get_beta_residues(self) -> List[int]:
        """Get residue numbers in beta strands."""
        return [r.resnum for r in self.residues.values() if r.ss == 'E']

    def get_helix_residues(self) -> List[int]:
        """Get residue numbers in helices (H, G, I)."""
        return [r.resnum for r in self.residues.values() if r.ss in 'HGI']


def run_dssp(structure_path: str, dssp_cmd: str = "mkdssp") -> Optional[str]:
    """
    Run DSSP on a structure file.

    Parameters:
    -----------
    structure_path : str
        Path to PDB or mmCIF file
    dssp_cmd : str
        DSSP command (mkdssp or dssp)

    Returns:
    --------
    str or None
        Path to DSSP output file, or None if DSSP failed
    """
    structure_path = Path(structure_path)

    # Create output path
    dssp_path = structure_path.with_suffix('.dssp')

    try:
        result = subprocess.run(
            [dssp_cmd, str(structure_path), str(dssp_path)],
            capture_output=True,
            text=True,
            timeout=60,
        )

        if result.returncode != 0:
            logger.warning(f"DSSP failed: {result.stderr}")
            return None

        if not dssp_path.exists():
            logger.warning("DSSP did not produce output file")
            return None

        return str(dssp_path)

    except FileNotFoundError:
        logger.warning(f"DSSP command '{dssp_cmd}' not found. Install with: conda install -c salilab dssp")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("DSSP timed out")
        return None
    except Exception as e:
        logger.warning(f"Error running DSSP: {e}")
        return None


def parse_dssp(dssp_path: str) -> Optional[DSSPData]:
    """
    Parse DSSP output file.

    Parameters:
    -----------
    dssp_path : str
        Path to DSSP output file

    Returns:
    --------
    DSSPData or None
        Parsed DSSP data, or None if parsing failed
    """
    residues = {}
    all_hbonds = []

    try:
        with open(dssp_path) as f:
            in_data = False

            for line in f:
                # Look for start of residue data
                if line.startswith('  #  RESIDUE'):
                    in_data = True
                    continue

                if not in_data:
                    continue

                # Skip if line is too short
                if len(line) < 75:
                    continue

                # Parse residue line
                try:
                    # Residue number (columns 1-5)
                    resnum_str = line[0:5].strip()
                    if not resnum_str or resnum_str == '!':
                        continue  # Chain break
                    resnum = int(resnum_str)

                    # Chain ID (column 12)
                    chain = line[11:12].strip() or 'A'

                    # Amino acid (column 14)
                    aa = line[13:14]
                    if aa == '!':
                        continue  # Chain break

                    # Secondary structure (column 17)
                    ss = line[16:17]
                    if ss == ' ':
                        ss = 'C'  # Coil

                    # Solvent accessibility (columns 35-38)
                    try:
                        acc = int(line[34:38].strip())
                    except ValueError:
                        acc = 0

                    # Parse hydrogen bonds
                    # Format: offset,energy pairs at specific column positions
                    # NH-->O bonds at columns 39-50 and 51-62
                    # O-->HN bonds at columns 63-74 (but these are redundant)

                    hbonds = []

                    # First NH-->O bond (columns 39-50)
                    hb1 = _parse_hbond_field(line[39:50], resnum, is_donor=True)
                    if hb1:
                        hbonds.append(hb1)
                        all_hbonds.append(hb1)

                    # Second NH-->O bond (columns 51-62)
                    hb2 = _parse_hbond_field(line[50:61], resnum, is_donor=True)
                    if hb2:
                        hbonds.append(hb2)
                        all_hbonds.append(hb2)

                    residues[resnum] = DSSPResidue(
                        resnum=resnum,
                        chain=chain,
                        aa=aa,
                        ss=ss,
                        acc=acc,
                        hbonds=hbonds,
                    )

                except (ValueError, IndexError) as e:
                    # Skip malformed lines
                    continue

        # Remove duplicate H-bonds (same pair might be listed twice)
        unique_hbonds = []
        seen = set()
        for hb in all_hbonds:
            key = (min(hb.donor_resnum, hb.acceptor_resnum),
                   max(hb.donor_resnum, hb.acceptor_resnum))
            if key not in seen:
                seen.add(key)
                unique_hbonds.append(hb)

        logger.info(f"Parsed DSSP: {len(residues)} residues, {len(unique_hbonds)} H-bonds "
                   f"({len([h for h in unique_hbonds if h.is_long_range])} long-range)")

        return DSSPData(residues=residues, hbonds=unique_hbonds)

    except Exception as e:
        logger.warning(f"Error parsing DSSP file: {e}")
        return None


def _parse_hbond_field(field: str, resnum: int, is_donor: bool) -> Optional[HydrogenBond]:
    """
    Parse a hydrogen bond field from DSSP.

    Format: "  -3,-0.3" means partner at resnum-3 with energy -0.3
    """
    try:
        field = field.strip()
        if not field or ',' not in field:
            return None

        parts = field.split(',')
        if len(parts) != 2:
            return None

        offset = int(parts[0].strip())
        energy = float(parts[1].strip())

        # Skip if no bond (offset 0 or very weak)
        if offset == 0 or energy > -0.1:
            return None

        partner_resnum = resnum + offset

        if is_donor:
            # This residue's NH donates to partner's O
            return HydrogenBond(
                donor_resnum=resnum,
                acceptor_resnum=partner_resnum,
                energy=energy,
            )
        else:
            # This residue's O accepts from partner's NH
            return HydrogenBond(
                donor_resnum=partner_resnum,
                acceptor_resnum=resnum,
                energy=energy,
            )

    except (ValueError, IndexError):
        return None


def get_dssp_for_structure(
    structure_path: str,
    cache_dir: Optional[Path] = None,
    dssp_cmd: str = "mkdssp",
) -> Optional[DSSPData]:
    """
    Get DSSP data for a structure, running DSSP if needed.

    DSSP output is cached in the cache_dir for future use.

    Parameters:
    -----------
    structure_path : str
        Path to PDB or mmCIF file
    cache_dir : Path, optional
        Directory to cache DSSP output
    dssp_cmd : str
        DSSP command to use

    Returns:
    --------
    DSSPData or None
    """
    structure_path = Path(structure_path)

    # Determine where to store/look for DSSP file
    if cache_dir:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        dssp_path = cache_dir / f"{structure_path.stem}.dssp"
    else:
        dssp_path = structure_path.with_suffix('.dssp')

    # Check for cached DSSP
    if dssp_path.exists():
        logger.info(f"Using cached DSSP: {dssp_path}")
        return parse_dssp(str(dssp_path))

    # Run DSSP and save to cache location
    logger.info(f"Running DSSP on {structure_path}")
    try:
        result = subprocess.run(
            [dssp_cmd, str(structure_path), str(dssp_path)],
            capture_output=True,
            text=True,
            timeout=60,
        )

        if result.returncode != 0:
            logger.warning(f"DSSP failed: {result.stderr}")
            return None

        if not dssp_path.exists():
            logger.warning("DSSP did not produce output file")
            return None

        logger.info(f"DSSP output cached to: {dssp_path}")
        return parse_dssp(str(dssp_path))

    except FileNotFoundError:
        logger.warning(f"DSSP command '{dssp_cmd}' not found. Install with: conda install -c salilab dssp")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("DSSP timed out")
        return None
    except Exception as e:
        logger.warning(f"Error running DSSP: {e}")
        return None
