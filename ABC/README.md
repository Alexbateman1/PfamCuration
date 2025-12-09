# ABC (Alex Bateman Chop) Domain Predictor

A geometric approach for predicting protein domains from AlphaFold structures using contact graph analysis and community detection.

## Hypothesis

Simple geometric principles (contact density, spatial compactness) can effectively identify protein domains without complex deep learning methods.

## Algorithm Overview

### 1. Build Contact Graph
- **Nodes**: Cα atoms (one per residue)
- **Edges**: Spatial proximity (residues within distance threshold)
- **Edge weights**: Combined function of distance, pLDDT confidence, and PAE

Weight function:
```
w_ij = (pLDDT_i × pLDDT_j / 10000) × exp(-d²_ij/σ²)
```

If PAE available:
```
w_ij *= exp(-PAE_ij / pae_scale)
```

### 2. Graph Clustering
- Community detection using Leiden or Louvain algorithm
- Resolution parameter controls granularity
- Clusters at multiple thresholds for stability analysis

### 3. Quality Assessment per Cluster
- **Radius of gyration** (compactness)
- **Contact density ratio**: internal / external contacts
- **Average pLDDT** (confidence)
- **Minimum size threshold** (~30 residues)

### 4. Refinement
- Merge over-split domains with high inter-cluster contact density
- Filter clusters below minimum size
- Optimize boundaries by reassigning edge residues

### 5. NDR Identification
- Residues not assigned to any cluster
- Low pLDDT regions (<70)
- Isolated residues with few contacts
- Terminal disordered regions

## Installation

### Requirements
```bash
pip install biopython networkx numpy scipy matplotlib

# For Leiden clustering (recommended)
pip install leidenalg igraph
```

### As part of PfamCuration
The ABC module is included in the PfamCuration repository.

## Usage

### Command Line

```bash
# Predict from UniProt accession
python -m ABC.cli predict --uniprot P12345

# Predict from local structure file
python -m ABC.cli predict --pdb structure.cif --pae pae.json

# Batch processing
python -m ABC.cli batch --input proteins.txt --output results/

# Custom parameters
python -m ABC.cli predict --uniprot P12345 \
    --distance-threshold 12 \
    --resolution 0.8 \
    --min-domain-size 25 \
    --ndr-cutoff 70
```

### Python API

```python
from ABC import ABCPredictor

# Initialize predictor
predictor = ABCPredictor(
    distance_threshold=10.0,   # Å
    min_domain_size=30,        # residues
    ndr_plddt_cutoff=70.0,     # pLDDT threshold for NDR
    clustering_method="leiden",
    resolution=1.0,
)

# Predict from UniProt
prediction = predictor.predict_from_uniprot("P12345")

# Or from local files
prediction = predictor.predict_from_file(
    "structure.cif",
    "pae.json",  # optional
)

# View results
print(prediction.summary())

# Generate visualization and ChimeraX commands
chimerax_commands = predictor.visualize(prediction, "output_prefix")
```

### Output Files

For each prediction, the following files are generated:
- `*.cxc` - ChimeraX commands for coloring by domain
- `*.pml` - PyMOL script for coloring by domain
- `*.png` - Domain architecture diagram
- `*.html` - Detailed HTML report

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `distance_threshold` | 10.0 Å | Maximum Cα-Cα distance for contact |
| `min_domain_size` | 30 | Minimum residues for valid domain |
| `ndr_plddt_cutoff` | 70.0 | pLDDT below which residues may be NDR |
| `clustering_method` | "leiden" | Clustering algorithm ("leiden" or "louvain") |
| `resolution` | 1.0 | Clustering resolution (higher = more clusters) |
| `sigma` | 8.0 | Gaussian decay parameter for edge weights |
| `use_pae` | True | Use PAE for edge weighting if available |

### Parameter Exploration

Test sensitivity at different thresholds:
```python
# Multi-threshold analysis
results = predictor.multi_threshold_analysis(
    "structure.cif",
    "pae.json",
    thresholds=[8.0, 10.0, 12.0, 15.0]
)

# Domains stable across thresholds are more reliable
for threshold, pred in results.items():
    print(f"{threshold}Å: {len(pred.domains)} domains")
```

## Output Format

### Domain Object
```python
@dataclass
class Domain:
    domain_id: int
    segments: List[Tuple[int, int]]  # [(start, end), ...]
    residue_indices: List[int]
    quality_metrics: Dict

# Example
domain.to_chopping_string()  # "24-111_266-345" (discontinuous)
domain.size                   # 167 residues
domain.is_discontinuous       # True
```

### Quality Metrics
```python
{
    "radius_of_gyration": 15.2,      # Å
    "contact_density_ratio": 3.5,    # internal/external
    "avg_plddt": 85.3,
    "min_plddt": 62.1,
    "quality_score": 78.5,           # 0-100 composite score
    "n_segments": 2,
    "sequence_coverage": 0.85,
}
```

### NDR Region
```python
@dataclass
class NDRRegion:
    start: int
    end: int
    residue_indices: List[int]
    avg_plddt: float
    reason: str  # "low_plddt", "isolated", "unassigned"
```

## Testing

```bash
# Run on test set
python ABC/test_abc.py --all-tests

# Test single protein
python ABC/test_abc.py --uniprot P12345

# Multi-threshold stability test
python ABC/test_abc.py --uniprot P12345 --multi-threshold

# Parameter sensitivity analysis
python ABC/test_abc.py --uniprot P12345 --sensitivity
```

## Comparison with Other Methods

ABC aims to complement existing domain prediction methods:

| Method | Approach | Strengths |
|--------|----------|-----------|
| **ABC** | Geometric/graph-based | Interpretable, fast, leverages pLDDT/PAE |
| Chainsaw | Sequence-based | Works without structure |
| Merizo | Structure-informed ML | High accuracy on benchmarks |
| UniDoc | Deep learning | Handles complex architectures |

## Architecture

```
ABC/
├── __init__.py          # Package exports
├── abc_predictor.py     # Main predictor class
├── contact_graph.py     # Contact graph builder
├── domain_quality.py    # Quality assessment
├── ndr_detector.py      # NDR identification
├── visualize.py         # Visualization outputs
├── cli.py               # Command line interface
├── test_abc.py          # Test suite
└── README.md            # This file
```

## Future Improvements

- [ ] Benchmark against CATH/SCOP domain boundaries
- [ ] Compare with Merizo/UniDoc/Chainsaw predictions
- [ ] Add secondary structure analysis (DSSP)
- [ ] Implement domain type classification
- [ ] Add support for multi-chain complexes
- [ ] Web interface integration

## Author

Alex Bateman (ABC - Alex Bateman Chop)

## License

Part of the PfamCuration toolkit.
