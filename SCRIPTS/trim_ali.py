#!/usr/bin/env python3
"""
Trim ragged alignment ends based on column occupancy.

This script analyzes alignment column occupancy (fraction of amino acids vs gaps)
and trims the N and C termini to produce cleaner alignments.

Usage:
    trim_ali.py <alignment_file> [options]

Examples:
    trim_ali.py SEED                           # Use default threshold (0.5)
    trim_ali.py SEED --threshold 0.6           # Use 60% occupancy threshold
    trim_ali.py SEED --method gradient         # Use max gradient method
    trim_ali.py SEED --report                  # Generate occupancy report
    trim_ali.py SEED --dry-run                 # Show what would be trimmed without output
"""

import argparse
import sys
import re
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass


@dataclass
class AlignedSequence:
    """Represents a sequence in the alignment."""
    full_id: str          # e.g., "A0A068TQG2.1/265-412"
    base_id: str          # e.g., "A0A068TQG2.1"
    orig_start: int       # Original start coordinate
    orig_end: int         # Original end coordinate
    sequence: str         # Aligned sequence with gaps


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Trim ragged alignment ends based on column occupancy.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Methods:
  threshold  - Trim to the first/last column where occupancy >= threshold (default)
  gradient   - Trim to the column with the maximum occupancy increase

Examples:
  %(prog)s SEED                        # Use default threshold of 0.5
  %(prog)s SEED -t 0.6                 # Use 60%% occupancy threshold
  %(prog)s SEED --method gradient      # Use gradient-based trimming
  %(prog)s SEED --report               # Show occupancy statistics
  %(prog)s SEED --dry-run --report     # Analyze without writing output
        """
    )
    parser.add_argument(
        'alignment_file',
        help='Path to the alignment file (e.g., SEED)'
    )
    parser.add_argument(
        '-t', '--threshold',
        type=float,
        default=0.5,
        help='Occupancy threshold for trimming (0.0-1.0, default: 0.5)'
    )
    parser.add_argument(
        '-m', '--method',
        choices=['threshold', 'gradient'],
        default='threshold',
        help='Trimming method (default: threshold)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file (default: stdout)'
    )
    parser.add_argument(
        '--report',
        action='store_true',
        help='Print occupancy report to stderr'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Analyze and report but do not output trimmed alignment'
    )
    parser.add_argument(
        '--min-occupancy',
        type=float,
        default=0.1,
        help='Minimum occupancy to consider (used with gradient method, default: 0.1)'
    )
    parser.add_argument(
        '--window',
        type=int,
        default=1,
        help='Smoothing window size for gradient calculation (default: 1)'
    )

    return parser.parse_args()


def parse_alignment(filename: str) -> List[AlignedSequence]:
    """
    Parse alignment file in Stockholm/Pfam format.

    Args:
        filename: Path to alignment file

    Returns:
        List of AlignedSequence objects
    """
    sequences = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if not line or line.startswith('#'):
                continue

            # Match lines like: A0A068TQG2.1/265-412 -------RKLALYMVVCTLITPFLVYKYLDCAPP...
            # Note: no $ anchor to allow trailing whitespace
            match = re.match(r'^(\S+)/(\d+)-(\d+)\s+(\S+)', line)
            if match:
                base_id = match.group(1)
                orig_start = int(match.group(2))
                orig_end = int(match.group(3))
                sequence = match.group(4)
                full_id = f"{base_id}/{orig_start}-{orig_end}"

                sequences.append(AlignedSequence(
                    full_id=full_id,
                    base_id=base_id,
                    orig_start=orig_start,
                    orig_end=orig_end,
                    sequence=sequence
                ))

    return sequences


def calculate_occupancy(sequences: List[AlignedSequence]) -> List[float]:
    """
    Calculate per-column occupancy (fraction of amino acids vs gaps).

    Args:
        sequences: List of aligned sequences

    Returns:
        List of occupancy values (0.0 to 1.0) for each column
    """
    if not sequences:
        return []

    alignment_length = len(sequences[0].sequence)
    num_sequences = len(sequences)
    occupancy = []

    for col in range(alignment_length):
        aa_count = 0
        for seq in sequences:
            char = seq.sequence[col] if col < len(seq.sequence) else '-'
            if char not in '-.':  # Count anything that's not a gap
                aa_count += 1
        occupancy.append(aa_count / num_sequences)

    return occupancy


def smooth_occupancy(occupancy: List[float], window: int) -> List[float]:
    """
    Apply simple moving average smoothing to occupancy values.

    Args:
        occupancy: Raw occupancy values
        window: Window size for smoothing

    Returns:
        Smoothed occupancy values
    """
    if window <= 1:
        return occupancy

    smoothed = []
    half_window = window // 2

    for i in range(len(occupancy)):
        start = max(0, i - half_window)
        end = min(len(occupancy), i + half_window + 1)
        smoothed.append(sum(occupancy[start:end]) / (end - start))

    return smoothed


def find_trim_position_threshold(
    occupancy: List[float],
    threshold: float,
    from_start: bool = True
) -> int:
    """
    Find the first position where occupancy meets or exceeds threshold.

    Args:
        occupancy: Per-column occupancy values
        threshold: Occupancy threshold (0.0-1.0)
        from_start: If True, search from start; if False, search from end

    Returns:
        Column index for trim position
    """
    if from_start:
        for i, occ in enumerate(occupancy):
            if occ >= threshold:
                return i
        return 0
    else:
        for i in range(len(occupancy) - 1, -1, -1):
            if occupancy[i] >= threshold:
                return i
        return len(occupancy) - 1


def find_trim_position_gradient(
    occupancy: List[float],
    min_occupancy: float = 0.1,
    from_start: bool = True,
    high_occupancy_stop: float = 0.9
) -> int:
    """
    Find the position with maximum occupancy increase (steepest gradient).

    The search stops once occupancy reaches high_occupancy_stop to avoid
    finding internal gaps.

    Args:
        occupancy: Per-column occupancy values
        min_occupancy: Minimum occupancy to consider as candidate
        from_start: If True, search from start; if False, search from end
        high_occupancy_stop: Stop searching once this occupancy is reached

    Returns:
        Column index for trim position
    """
    if len(occupancy) < 2:
        return 0 if from_start else len(occupancy) - 1

    # Calculate gradients (differences between adjacent columns)
    gradients = []
    for i in range(len(occupancy) - 1):
        if from_start:
            grad = occupancy[i + 1] - occupancy[i]
        else:
            grad = occupancy[i] - occupancy[i + 1]
        gradients.append(grad)

    # For N-terminus: search from start until we hit high occupancy
    # For C-terminus: search from end until we hit high occupancy
    if from_start:
        search_end = 0
        for i, occ in enumerate(occupancy):
            if occ >= high_occupancy_stop:
                search_end = i
                break
        search_end = max(search_end, 1)  # At least search first column
        search_range = range(search_end)
    else:
        search_start = len(occupancy) - 1
        for i in range(len(occupancy) - 1, -1, -1):
            if occupancy[i] >= high_occupancy_stop:
                search_start = i
                break
        search_start = min(search_start, len(gradients) - 1)
        search_range = range(search_start, len(gradients))

    # Find position with maximum gradient where occupancy exceeds minimum
    best_pos = 0 if from_start else len(occupancy) - 1
    best_grad = -float('inf')

    for i in search_range:
        if i >= len(gradients):
            continue
        if gradients[i] > best_grad:
            # For N-terminus, check that the position after has sufficient occupancy
            # For C-terminus, check the position before
            check_pos = i + 1 if from_start else i
            if occupancy[check_pos] >= min_occupancy:
                best_grad = gradients[i]
                best_pos = i + 1 if from_start else i

    return best_pos


def build_position_mappings(seq: AlignedSequence) -> Tuple[List[Optional[int]], List[Optional[int]]]:
    """
    Build bidirectional position mappings between sequence and alignment positions.

    Args:
        seq: Aligned sequence

    Returns:
        Tuple of (seq_to_aln, aln_to_seq) mappings
    """
    seq_to_aln = []  # Index by (seq_pos - orig_start), value is alignment position
    aln_to_seq = []  # Index by alignment position, value is sequence position

    seq_pos = seq.orig_start

    for aln_pos, char in enumerate(seq.sequence):
        if char not in '-.':  # It's a residue
            while len(seq_to_aln) < (seq_pos - seq.orig_start):
                seq_to_aln.append(None)
            seq_to_aln.append(aln_pos)
            aln_to_seq.append(seq_pos)
            seq_pos += 1
        else:
            aln_to_seq.append(None)

    return seq_to_aln, aln_to_seq


def trim_alignment(
    sequences: List[AlignedSequence],
    trim_start: int,
    trim_end: int
) -> List[Tuple[str, str]]:
    """
    Trim alignment to specified column range and update coordinates.

    Args:
        sequences: List of aligned sequences
        trim_start: Start column (0-indexed, inclusive)
        trim_end: End column (0-indexed, inclusive)

    Returns:
        List of (new_header, trimmed_sequence) tuples
    """
    results = []

    for seq in sequences:
        # Build position mappings for this sequence
        seq_to_aln, aln_to_seq = build_position_mappings(seq)

        # Extract the trimmed region
        trimmed_seq = seq.sequence[trim_start:trim_end + 1]

        # Find new start position (first non-gap in trimmed region)
        new_start = None
        for i in range(trim_start, trim_end + 1):
            if i < len(aln_to_seq) and aln_to_seq[i] is not None:
                new_start = aln_to_seq[i]
                break

        # Find new end position (last non-gap in trimmed region)
        new_end = None
        for i in range(trim_end, trim_start - 1, -1):
            if i < len(aln_to_seq) and aln_to_seq[i] is not None:
                new_end = aln_to_seq[i]
                break

        # Handle all-gap sequences
        if new_start is None or new_end is None:
            new_start = seq.orig_start
            new_end = seq.orig_start - 1  # Empty range indication

        new_header = f"{seq.base_id}/{new_start}-{new_end}"
        results.append((new_header, trimmed_seq))

    return results


def print_occupancy_report(
    occupancy: List[float],
    trim_start: int,
    trim_end: int,
    method: str,
    threshold: float,
    file=sys.stderr
):
    """
    Print a report showing occupancy statistics and trim positions.
    """
    print("\n" + "=" * 70, file=file)
    print("ALIGNMENT OCCUPANCY REPORT", file=file)
    print("=" * 70, file=file)

    print(f"\nAlignment length: {len(occupancy)} columns", file=file)
    print(f"Trimming method: {method}", file=file)
    if method == 'threshold':
        print(f"Threshold: {threshold:.1%}", file=file)

    print(f"\nTrim positions:", file=file)
    print(f"  N-terminus: column {trim_start} (trimming {trim_start} columns)", file=file)
    print(f"  C-terminus: column {trim_end} (trimming {len(occupancy) - trim_end - 1} columns)", file=file)
    print(f"  New length: {trim_end - trim_start + 1} columns", file=file)

    # Show occupancy at key positions
    print(f"\nOccupancy at trim positions:", file=file)
    if trim_start > 0:
        print(f"  Before N-trim (col {trim_start-1}): {occupancy[trim_start-1]:.1%}", file=file)
    print(f"  At N-trim (col {trim_start}): {occupancy[trim_start]:.1%}", file=file)
    print(f"  At C-trim (col {trim_end}): {occupancy[trim_end]:.1%}", file=file)
    if trim_end < len(occupancy) - 1:
        print(f"  After C-trim (col {trim_end+1}): {occupancy[trim_end+1]:.1%}", file=file)

    # Summary statistics
    trimmed_occupancy = occupancy[trim_start:trim_end + 1]
    if trimmed_occupancy:
        print(f"\nTrimmed alignment statistics:", file=file)
        print(f"  Mean occupancy: {sum(trimmed_occupancy)/len(trimmed_occupancy):.1%}", file=file)
        print(f"  Min occupancy: {min(trimmed_occupancy):.1%}", file=file)
        print(f"  Max occupancy: {max(trimmed_occupancy):.1%}", file=file)

    # ASCII visualization of occupancy
    print(f"\nOccupancy profile (N -> C):", file=file)
    print("  0%       50%      100%", file=file)
    print("  |--------|--------|", file=file)

    # Show first and last 20 columns with markers
    n_show = 20

    def make_bar(occ, is_trimmed=False):
        bar_len = int(occ * 20)
        bar = '#' * bar_len + '.' * (20 - bar_len)
        marker = ' ' if is_trimmed else '*'
        return f"{marker} {bar} {occ:.0%}"

    print("\n  First columns:", file=file)
    for i in range(min(n_show, len(occupancy))):
        is_trimmed = i >= trim_start
        print(f"  {i:3d}: {make_bar(occupancy[i], is_trimmed)}", file=file)

    if len(occupancy) > 2 * n_show:
        print(f"  ... ({len(occupancy) - 2*n_show} columns) ...", file=file)

    print("\n  Last columns:", file=file)
    start_idx = max(n_show, len(occupancy) - n_show)
    for i in range(start_idx, len(occupancy)):
        is_trimmed = i <= trim_end
        print(f"  {i:3d}: {make_bar(occupancy[i], is_trimmed)}", file=file)

    print("\n  Legend: * = trimmed away, # = occupancy", file=file)
    print("=" * 70 + "\n", file=file)


def main():
    args = parse_arguments()

    # Parse the alignment
    try:
        sequences = parse_alignment(args.alignment_file)
    except FileNotFoundError:
        print(f"Error: File not found: {args.alignment_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing alignment: {e}", file=sys.stderr)
        sys.exit(1)

    if not sequences:
        print("Error: No sequences found in alignment", file=sys.stderr)
        sys.exit(1)

    # Calculate occupancy
    occupancy = calculate_occupancy(sequences)

    if args.method == 'gradient' and args.window > 1:
        occupancy_for_trim = smooth_occupancy(occupancy, args.window)
    else:
        occupancy_for_trim = occupancy

    # Find trim positions
    if args.method == 'threshold':
        trim_start = find_trim_position_threshold(occupancy_for_trim, args.threshold, from_start=True)
        trim_end = find_trim_position_threshold(occupancy_for_trim, args.threshold, from_start=False)
    else:  # gradient
        trim_start = find_trim_position_gradient(
            occupancy_for_trim,
            min_occupancy=args.min_occupancy,
            from_start=True
        )
        trim_end = find_trim_position_gradient(
            occupancy_for_trim,
            min_occupancy=args.min_occupancy,
            from_start=False
        )

    # Validate trim positions
    if trim_start >= trim_end:
        print(f"Warning: Invalid trim positions (start={trim_start}, end={trim_end}). ",
              file=sys.stderr)
        print("Try adjusting threshold or method.", file=sys.stderr)
        trim_start = 0
        trim_end = len(occupancy) - 1

    # Print report if requested
    if args.report or args.dry_run:
        print_occupancy_report(
            occupancy, trim_start, trim_end,
            args.method, args.threshold
        )

    # Output trimmed alignment
    if not args.dry_run:
        trimmed = trim_alignment(sequences, trim_start, trim_end)

        # Determine output destination
        if args.output:
            out_file = open(args.output, 'w')
        else:
            out_file = sys.stdout

        try:
            for header, seq in trimmed:
                print(f"{header:25s} {seq}", file=out_file)
        finally:
            if args.output:
                out_file.close()


if __name__ == '__main__':
    main()
