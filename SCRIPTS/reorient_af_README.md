# reorient_af.py

Reorient AlphaFold structures to minimize visual overlap of TED domains.

## Overview

This script automatically:
1. Downloads AlphaFold structural models from AlphaFold DB (in CIF format)
2. Fetches domain definitions from the TED database
3. Calculates an optimal rotation to minimize visual overlap (occlusion) of domains
4. Outputs a reoriented structure file with domain annotations
5. Creates a 3D visualization with TED domains colored according to the official TED color scheme

## Installation

### Dependencies

```bash
pip install -r reorient_af_requirements.txt
```

Required packages:
- biopython >= 1.79
- numpy >= 1.20.0
- matplotlib >= 3.3.0

## Usage

### Basic usage:

```bash
./reorient_af.py <UniProt_accession>
```

### Examples:

```bash
# Simple usage with default output files
./reorient_af.py P12345

# Specify custom output file names
./reorient_af.py P12345 --output-cif my_structure.cif --output-png my_viz.png

# With version number
./reorient_af.py P12345.2

# Use more samples for better optimization (slower but more thorough)
./reorient_af.py P12345 --n-samples 5000
```

### Output Files

The script generates two files:

1. **ted_reoriented.cif** (default) - The reoriented structure in CIF format with TED domain annotations in the header
2. **ted3d.png** (default) - A 3D visualization of the structure with colored TED domains

## Algorithm

The script uses **rotation sampling with 2D overlap minimization** to find the optimal viewing orientation:

1. Extract C-alpha atoms for each TED domain
2. Sample 1000 viewing angles using Fibonacci sphere sampling for uniform coverage
3. For each angle:
   - Rotate the structure
   - Project C-alphas to 2D (XY plane)
   - Calculate bounding box overlap between domain pairs
4. Select the rotation with minimal overlap

The algorithm uses only C-alpha atoms for calculations (much faster), but rotates all atoms for the output structure. The visualization also shows only the C-alpha backbone trace.

## TED Domain Color Scheme

The script uses the official TED database color scheme:

| Domain | Color   | Hex Code |
|--------|---------|----------|
| 01     | Blue    | #4A79A7  |
| 02     | Orange  | #F28E2C  |
| 03     | Red     | #E15759  |
| 04     | Teal    | #76B7B2  |
| 05     | Green   | #59A14F  |
| 06     | Yellow  | #EDC949  |
| 07     | Purple  | #AF7AA1  |
| 08     | Pink    | #FF9DA7  |
| 09     | Brown   | #9C755F  |
| 10     | Grey    | #BAB0AB  |

Colors cycle if there are more than 10 domains.

## Notes

- The script requires internet access to fetch data from AlphaFold DB and TED
- If no TED domains are found for the protein, the script will exit with an error
- The rotation is optimized for viewing, not for any biological interpretation
- Domain annotations are included as comments in the output CIF file

## API Endpoints Used

- **AlphaFold DB**: `https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.cif`
- **TED**: `https://ted.cathdb.info/api/v1/uniprot/summary/{accession}`
