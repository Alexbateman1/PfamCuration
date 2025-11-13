# Pfam Curation Codebase Exploration Summary

## Overview
The PfamCuration repository contains a comprehensive toolkit for Pfam database curation, focusing on building and refining protein family profiles from UniProt data.

---

## 1. EXISTING PIPELINE INFRASTRUCTURE

### Main Pipeline: pfam_curation_pipeline.py (22 KB, ~620 lines)
**Location:** `/home/user/PfamCuration/SCRIPTS/pfam_curation_pipeline.py`

**Purpose:** Automated end-to-end curation pipeline that orchestrates multiple curation stages.

**Key Features:**
- **Input Methods:** UniProt proteome IDs or accession lists
- **Progress Tracking:** Hidden marker files (.ted_build_complete, .iteration_complete, etc.)
- **Resume/Reset:** Supports resuming from last completed step or resetting progress
- **SLURM Integration:** Monitors job completion via `squeue`

**Pipeline Stages:**
1. **TED Build** - Creates domains from UniProt accessions using TED API
2. **Initial Triage** - Runs triage.pl to analyze family quality
3. **Iteration Cycles** - Refines alignments until no new improvement directories
4. **Annotation Scripts** - Adds references, PDB links, SwissProt data, etc.

---

## 2. CORE ALIGNMENT & FAMILY BUILDING

### iterate.py (44 KB, ~1000+ lines)
**Location:** `/home/user/PfamCuration/SCRIPTS/iterate.py`

**Purpose:** Python port of iterate.pl for iterative family refinement.

**Key Algorithms & Functions:**

#### SEED Alignment Construction Methods:

1. **SEEDify Methodology** (Advanced):
   - Uses best-scoring sequence as template
   - Extracts SwissProt sequences
   - Applies `extend.pl` for alignment extension
   - Uses `split_align.pl` to trim to template coordinates
   - Uses `belvu` for partial sequence removal

2. **Traditional SEED Building** (Fallback):
   - Progressive non-redundancy reduction with `belvu`
   - Starts at 80% identity, reduces by 10% per iteration
   - Uses `create_alignment.pl` with MAFFT realignment

#### Key External Tools Used:
- **`pfco`** - Pfam checkout/version control
- **`belvu`** - Alignment visualization & processing (non-redundancy)
- **`extend.pl`** - Extends alignments while preserving gaps
- **`split_align.pl`** - Trims alignments to region coordinates
- **`create_alignment.pl`** - Alignment creation with MAFFT
- **`pfbuild`** - Builds HMM profiles
- **`pqc-overlap-rdb`** - Checks overlap with existing Pfam families
- **`swissprot.pl`** - Extracts SwissProt sequences

#### Quality Control Features:
- DNI (Do Not Iterate) flag detection
- HT (High Threshold) flag detection
- Partial sequence removal for non-repeat families
- Repeat family special handling
- Alignment hash tracking for benchmarking

---

## 3. TED DOMAIN PROCESSING

### ted_build.py (150+ lines)
**Location:** `/home/user/PfamCuration/SCRIPTS/ted_build.py`

**Purpose:** Creates curation directories from TED (Topology Extracted Domains) domains.

**Key Features:**
- Fetches TED domains from `https://ted.cathdb.info/api/v1/`
- Downloads PDB files for TED domains
- Checks overlap with existing Pfam domains using InterPro API
- Skips domains with >50% overlap
- FoldSeek database search options

**External API Integration:**
- TED CATHDB API
- InterPro API for Pfam domain lookup

---

## 4. ANNOTATION & REFERENCE SCRIPTS

### add_pdb_ref.py
**Purpose:** Identifies PDB structures matching family sequences
- Uses SIFTS (Structure Integration with Function, Taxonomy and Sequences) mapping
- Matches UniProt accessions to PDB codes
- Adds structure papers to DESC files

### add_abb_ref.py
**Purpose:** Queries UniProt additional bibliography
- Searches UniProt bibliography file for accessions in ALIGN
- Retrieves supplementary references
- Integrates with DESC file updates

### query_paperblast.py
**Purpose:** Searches PaperBLAST database using HMM profiles
- Uses `hmmsearch` to query profiles against database
- Parses hmmsearch tblout format
- Queries SQLite database for gene/paper associations
- Returns results sorted by E-value

### triage_helper.py
**Purpose:** Interactive triage directory analysis
- Counts sequences in ALIGN/SEED files
- Identifies best family versions
- Filters by SwissProt content
- Detects overlaps with existing families

### clan_network_viz.py (46 KB)
**Purpose:** Visualizes clan-level family relationships
- **MySQL Database Access:** Queries pfam_live database
- **FoldSeek Integration:** Processes structure search results
- **Interactive HTML Export:** Creates network visualizations with vis.js
- **Features:**
  - Node type classification (Family, Domain, Repeat)
  - E-value-based edge coloring
  - Family information display
  - Interactive network exploration

---

## 5. DATA STRUCTURES & FILE FORMATS

### Directory Organization
```
CURATION_DIRECTORY/
├── ALIGN              # Stockholm format alignment file
├── SEED               # Manually curated seed alignment
├── DESC               # Pfam family metadata
├── sp                 # SwissProt sequences
├── scores             # HMM search scores
├── overlap            # Overlap with existing families
├── triage             # Quality metrics file
└── Iterate/           # Iteration subdirectories
    └── [same structure]
```

### Key Files:
- **ALIGN**: Stockholm format multiple sequence alignment
- **SEED**: High-quality reference alignment
- **DESC**: Family description file with metadata
- **scores**: HMM-search output with sequence scores
- **triage**: Tab-separated metrics (sequences, overlaps, quality scores)

### Triage File Format (5+ columns):
```
directory_name  total_seqs  overlaps  non_overlap_seqs  quality_fraction  [additional_info...]
```

---

## 6. EXTERNAL TOOL ECOSYSTEM

### Pfam-Specific Tools:
- **pfco** - Pfam version control checkout
- **pfbuild** - Profile HMM builder
- **pqc-overlap-rdb** - Overlap detection against RDB
- **belvu** - Stockholm alignment viewer/processor
- **extend.pl** - Alignment extension preserving annotation
- **split_align.pl** - Region-specific alignment extraction
- **create_alignment.pl** - MAFFT-based alignment creation
- **swissprot.pl** - Swiss-Prot sequence extraction
- **add_author.pl** - Author metadata addition
- **ted_ali.pl** - TED annotation addition
- **species_summary.pl** - Species distribution analysis
- **add_swiss_ref.pl** - Swiss-Prot reference linking

### Bioinformatics Tools:
- **MAFFT** - Multiple sequence alignment (used via create_alignment.pl)
- **HMMER** (implied) - HMM operations
- **FoldSeek** - Structure similarity searching
- **BLAST** (via PaperBLAST database)

### Visualization:
- **vis.js** - Network visualization (in clan_network_viz.py)
- **Matplotlib** (implied from pandas usage)

---

## 7. API & DATABASE INTEGRATION

### Remote APIs:
1. **UniProt REST API**
   - Proteome retrieval: `https://rest.uniprot.org/uniprotkb/stream`
   - Accession lookup
   - 5-minute timeout for large proteomes

2. **TED CATHDB API**
   - Domain information: `/api/v1/uniprot/summary/`
   - PDB file download: `/api/v1/files/{ted_id}.pdb`

3. **InterPro API**
   - Pfam domain lookup: `/interpro/api/entry/interpro/protein/uniprot/`

### Local Databases:
- **Pfam MySQL Database** (pfam_live)
  - Clan information
  - Family metadata
  - Sequence information

- **PaperBLAST SQLite Database**
  - GenePaper table (text-mined papers)
  - Gene information
  - Snippet mapping

- **SIFTS Mapping File**
  - UniProt to PDB chain mapping
  - Located at: `/homes/agb/Scripts/pdb_chain_uniprot.csv`

- **UniProt Bibliography File**
  - Supplementary references
  - Accession-based lookup

---

## 8. SVN INTERACTION (MINIMAL/LEGACY)

### Current Status:
- **No active SVN integration found** in the codebase
- `.svn` directories are explicitly skipped in triage_merge.py:
  ```python
  skip_dirs = {'DONE', 'IGNORE', 'OVERLAP', '.git', '.svn'}
  ```
- Git is the current version control system

### Implications for New Pipeline:
- SVN integration would need to be newly implemented
- No existing SVN client calls or authentication mechanisms present
- Candidate for new feature addition

---

## 9. HMM & ALIGNMENT PROCESSING CAPABILITIES

### HMM Operations:
- **pfbuild** constructs HMM models from SEED alignments
- **extend.pl** handles HMM-based alignment extension
- **hmmsearch** (via query_paperblast.py) searches profiles

### Alignment Quality Assessment:
```python
is_flush(alignment_file)  # Fraction of sequences flush to alignment ends
num_seq(filename)         # Count sequences
get_alignment_hash()      # MD5 hash for change detection
```

### SEED Building Pipeline:
1. Filter non-redundant sequences (belvu)
2. Remove partial sequences
3. Extend with template guidance
4. Realign with MAFFT
5. Mark with .goodseed flag

---

## 10. CONFIGURATION & HARDCODED PATHS

### Known Hardcoded Paths:
- Iteration archive: `/homes/agb/Curation/ITERATION_ARCHIVE`
- Benchmark directory: `/nfs/production/agb/pfam/users/agb/ALI_BENCHMARK`
- iterate.pl: `/homes/agb/Scripts/iterate_inline.pl`
- SIFTS file: `/homes/agb/Scripts/pdb_chain_uniprot.csv`
- Clan script: `/nfs/production/agb/interpro/users/typhaine/interpro-pfam-curation-tools/...`

### MySQL Connection:
- Config file: `~/.my.cnf`
- Database: `pfam_live`
- Default host: localhost

---

## 11. WHAT'S MISSING / NOT PRESENT

### Critical Gaps for SEED Extraction Pipeline:
1. **SVN Repository Access** - Would need to implement:
   - SVN client integration (subprocess calls to `svn` commands)
   - Authentication mechanism
   - Repository path configuration
   - Caching of extracted alignments/HMMs

2. **Bulk HHblits Execution** - Not currently implemented:
   - No HHblits command integration found
   - No all-against-all comparison framework
   - No result aggregation/comparison logic

3. **Configuration System** - Currently uses:
   - Command-line arguments
   - Hardcoded paths (problematic for portability)
   - Individual script configurations
   - Would benefit from YAML/JSON config file

4. **Parallel Execution Framework** - Currently uses:
   - SLURM for pfbuild jobs
   - No integrated job submission for general scripts
   - Could be enhanced for HHblits batch processing

5. **Result Aggregation & Analysis** - Not present:
   - No summary statistics compilation
   - No comparison result ranking
   - No visualization of all-against-all results

---

## 12. RELEVANT FOR NEW PIPELINE

### Reusable Components:
1. ✓ **UniProt API integration** (pfam_curation_pipeline.py)
   - Can reuse proteome fetching logic
   - Accession handling patterns

2. ✓ **SLURM job management** (pfam_curation_pipeline.py)
   - Job monitoring with `squeue`
   - Timeout handling
   - Job completion checking

3. ✓ **File I/O patterns**
   - Stockholm format handling
   - Tab-separated triage files
   - Progress tracking with hidden marker files

4. ✓ **Error handling & logging**
   - Established logging patterns
   - Subprocess error handling
   - Fallback mechanisms

5. ✓ **Alignment processing utilities** (iterate.py)
   - Sequence counting
   - Quality metrics
   - Template-based trimming logic

### Tools to Integrate With:
- **belvu** - For alignment viewing/filtering
- **MAFFT** (via create_alignment.pl) - For alignment generation
- **pfbuild** - For HMM profile creation
- **hmmsearch** - For profile searching

### Architecture to Consider:
- Directory structure per SEED/HMM
- Progress tracking with marker files
- Triage-style quality metrics file
- SLURM-based parallel execution
- MySQL database integration for Pfam metadata

---

## 13. SUMMARY TABLE OF SCRIPTS

| Script | Lines | Purpose | Key Dependencies |
|--------|-------|---------|-----------------|
| pfam_curation_pipeline.py | 620 | Orchestrate full curation pipeline | requests, SLURM (squeue) |
| iterate.py | 1000+ | Iterative family refinement | subprocess, belvu, pfbuild, pfco |
| ted_build.py | 150+ | TED domain discovery | requests, json, InterPro API |
| add_pdb_ref.py | 100+ | PDB structure linking | SIFTS file, urllib |
| add_abb_ref.py | 80+ | Bibliography integration | UniProt bibl file, subprocess |
| query_paperblast.py | 200+ | Paper discovery via HMM | sqlite3, subprocess (hmmsearch) |
| triage_helper.py | 150+ | Quality assessment | subprocess, file I/O |
| clan_network_viz.py | 1500+ | Network visualization | mysql.connector, pandas, vis.js |
| search_pfam_structure.py | 200+ | AlphaFold structure search | mysql.connector, FoldSeek |
| triage_merge.py | 300+ | Merge family identification | subprocess |

---

## 14. KEY RECOMMENDATIONS FOR HHBLITS PIPELINE

1. **SVN Access Layer:**
   - Create modular SVN wrapper class
   - Implement caching mechanism
   - Handle authentication via config

2. **SEED/HMM Extraction:**
   - Reuse pfam_curation_pipeline.py structure
   - Extend with SVN checkout step
   - Add SEED/HMM extraction methods

3. **HHblits Execution:**
   - Leverage SLURM job submission framework
   - Batch query generation
   - Progress tracking

4. **Result Management:**
   - Triage-style metrics file for comparisons
   - Database storage option
   - Visualization generation

5. **Configuration:**
   - Create config.yaml for paths, repositories
   - Support environment variable overrides
   - Centralize hardcoded paths

