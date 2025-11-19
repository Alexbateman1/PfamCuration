#!/usr/bin/env python3
"""
Simple test of HHsearch with direct HHM comparison (no database).
"""

import sys
import subprocess
from pathlib import Path

def test_family(family_id, seed_dir, work_dir):
    """Test conversion and hhmake only."""

    work_dir = Path(work_dir)
    work_dir.mkdir(exist_ok=True)

    seed_file = Path(seed_dir) / f"{family_id}_SEED"
    if not seed_file.exists():
        print(f"ERROR: SEED file not found: {seed_file}")
        return False

    print(f"Testing {family_id}...")

    # Step 1: Convert SEED to aligned FASTA
    a3m_file = work_dir / f"{family_id}.a3m"
    print(f"\n1. Converting SEED to aligned FASTA...")

    sequences = {}
    with open(seed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(None, 1)
            if len(parts) >= 2:
                seq_id = parts[0]
                seq = parts[1].replace('.', '-')
                sequences[seq_id] = seq

    with open(a3m_file, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")

    print(f"   ✓ Created: {a3m_file}")
    print(f"   ✓ Sequences: {len(sequences)}")

    # Check all sequences have same length
    lengths = set(len(seq) for seq in sequences.values())
    if len(lengths) > 1:
        print(f"   ✗ ERROR: Sequences have different lengths: {lengths}")
        return False
    else:
        print(f"   ✓ All sequences same length: {list(lengths)[0]}")

    # Step 2: Build HHM profile with hhmake
    hhm_file = work_dir / f"{family_id}.hhm"
    print(f"\n2. Building HHM profile with hhmake -M first...")

    cmd = f"hhmake -i {a3m_file} -o {hhm_file} -M first"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"   ✗ ERROR: hhmake failed")
        print(f"   STDERR: {result.stderr}")
        return False

    print(f"   ✓ Created: {hhm_file}")
    print(f"   ✓ Size: {hhm_file.stat().st_size} bytes")

    # Step 3: Test HHsearch with direct HHM comparison (self vs self)
    hhr_file = work_dir / f"{family_id}_vs_self.hhr"
    print(f"\n3. Running HHsearch (self-comparison test)...")

    # Use -t for target HHM (no database needed)
    cmd = f"hhsearch -i {hhm_file} -t {hhm_file} -o {hhr_file} -e 1.0"
    print(f"   Command: {cmd}")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)

    if result.returncode != 0:
        print(f"   ✗ ERROR: hhsearch failed")
        print(f"   STDERR: {result.stderr}")
        return False

    print(f"   ✓ Created: {hhr_file}")
    print(f"   ✓ Size: {hhr_file.stat().st_size} bytes")

    # Step 4: Parse results
    print(f"\n4. Parsing HHR results...")

    with open(hhr_file, 'r') as f:
        content = f.read()

    if family_id in content:
        print(f"   ✓ Found {family_id} in results (self-hit expected)")
    else:
        print(f"   ✗ WARNING: {family_id} not found in results")

    # Show first 30 lines
    lines = content.split('\n')
    print(f"\n5. HHR file preview (first 30 lines):")
    for i, line in enumerate(lines[:30], 1):
        print(f"   {i:3d}: {line}")

    print(f"\n{'='*60}")
    print(f"✅ SUCCESS! All steps completed for {family_id}")
    print(f"{'='*60}")
    print(f"\nFiles created in {work_dir}:")
    print(f"  - {a3m_file.name} (aligned FASTA)")
    print(f"  - {hhm_file.name} (HHM profile)")
    print(f"  - {hhr_file.name} (HHsearch results)")

    return True


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 test_simple.py <family_id> <seed_dir> <work_dir>")
        print("")
        print("Example:")
        print("  python3 test_simple.py PF00001 \\")
        print("    /nfs/production/agb/pfam/users/agb/HHsearch/DATA/SEED \\")
        print("    /tmp/test_hhsearch")
        sys.exit(1)

    family_id = sys.argv[1]
    seed_dir = sys.argv[2]
    work_dir = sys.argv[3]

    success = test_family(family_id, seed_dir, work_dir)
    sys.exit(0 if success else 1)
