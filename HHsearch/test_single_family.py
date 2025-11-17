#!/usr/bin/env python3
"""
Test HHsearch pipeline with a single family.
Quick sanity check before running the full pipeline.
"""

import sys
import subprocess
from pathlib import Path
import tempfile
import shutil

def test_family(family_id, seed_dir, work_dir):
    """Test the complete workflow for one family."""

    work_dir = Path(work_dir)
    work_dir.mkdir(exist_ok=True)

    seed_file = Path(seed_dir) / f"{family_id}_SEED"
    if not seed_file.exists():
        print(f"ERROR: SEED file not found: {seed_file}")
        return False

    print(f"Testing {family_id}...")
    print(f"SEED file: {seed_file}")

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

    print(f"   Created: {a3m_file}")
    print(f"   Sequences: {len(sequences)}")

    # Check all sequences have same length
    lengths = set(len(seq) for seq in sequences.values())
    if len(lengths) > 1:
        print(f"   WARNING: Sequences have different lengths: {lengths}")
    else:
        print(f"   All sequences same length: {list(lengths)[0]}")

    # Step 2: Build HHM profile with hhmake
    hhm_file = work_dir / f"{family_id}.hhm"
    print(f"\n2. Building HHM profile with hhmake -M first...")

    cmd = f"hhmake -i {a3m_file} -o {hhm_file} -M first"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"   ERROR: hhmake failed")
        print(f"   STDERR: {result.stderr}")
        return False

    print(f"   Created: {hhm_file}")
    print(f"   Size: {hhm_file.stat().st_size} bytes")

    # Step 3: Create a tiny database (just this one HHM)
    db_dir = work_dir / "test_db"
    db_dir.mkdir(exist_ok=True)

    print(f"\n3. Creating test database...")

    # Copy HHM to database directory
    shutil.copy(hhm_file, db_dir)

    # Build ffindex database
    cmd = f"ffindex_build -s {db_dir}/test_db.ffdata {db_dir}/test_db.ffindex {db_dir}/"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"   ERROR: ffindex_build failed")
        print(f"   STDERR: {result.stderr}")
        return False

    # Create symlinks for cs219 and a3m
    for suffix in ['cs219', 'a3m']:
        for ext in ['ffdata', 'ffindex']:
            target = db_dir / f"test_db.{ext}"
            link = db_dir / f"test_db_{suffix}.{ext}"
            if not link.exists() and target.exists():
                link.symlink_to(target.name)

    print(f"   Database: {db_dir}/test_db")

    # Step 4: Run HHsearch
    hhr_file = work_dir / f"{family_id}.hhr"
    print(f"\n4. Running HHsearch...")

    cmd = f"hhsearch -i {hhm_file} -d {db_dir}/test_db -o {hhr_file} -cpu 1 -e 1.0"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"   ERROR: hhsearch failed")
        print(f"   STDERR: {result.stderr}")
        return False

    print(f"   Created: {hhr_file}")
    print(f"   Size: {hhr_file.stat().st_size} bytes")

    # Step 5: Parse results
    print(f"\n5. Parsing HHR results...")

    with open(hhr_file, 'r') as f:
        lines = f.readlines()

    # Find hits
    in_table = False
    hit_count = 0
    for line in lines:
        if line.startswith(' No Hit'):
            in_table = True
            continue
        if in_table and line.strip() and not line.startswith('No ') and not line.startswith('Done'):
            hit_count += 1

    print(f"   Hits found: {hit_count}")

    # Show first few lines of HHR
    print(f"\n6. HHR file preview (first 20 lines):")
    for i, line in enumerate(lines[:20], 1):
        print(f"   {i:3d}: {line.rstrip()}")

    print(f"\n✅ SUCCESS! Test completed for {family_id}")
    return True


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 test_single_family.py <family_id> <seed_dir> <work_dir>")
        print("")
        print("Example:")
        print("  python3 test_single_family.py PF00001 \\")
        print("    /nfs/production/agb/pfam/users/agb/HHsearch/DATA/SEED \\")
        print("    /tmp/hhsearch_test")
        sys.exit(1)

    family_id = sys.argv[1]
    seed_dir = sys.argv[2]
    work_dir = sys.argv[3]

    success = test_family(family_id, seed_dir, work_dir)
    sys.exit(0 if success else 1)
