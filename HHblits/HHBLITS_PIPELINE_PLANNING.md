# HHBLITS Pipeline Planning Guide
## Quick Reference for SEED Alignment & HMM Extraction with All-Against-All Comparisons

### Key Findings from Codebase Analysis

---

## 1. IMMEDIATE GAPS TO ADDRESS

### SVN Integration (NEW REQUIREMENT)
**Current Status:** No SVN integration exists
**What's Needed:**
- SVN client wrapper class using `subprocess`
- Authentication configuration (via `~/.svn` credentials or command-line)
- Repository URL/path configuration
- Caching strategy for extracted files (to avoid redundant SVN calls)

**Example Structure:**
```python
class SVNManager:
    def __init__(self, repo_url, cache_dir):
        self.repo_url = repo_url
        self.cache_dir = cache_dir
    
    def checkout_family(self, family_id, output_dir):
        # svn checkout <repo_url>/families/<family_id> <output_dir>
        pass
    
    def get_seed_alignment(self, family_id):
        # Extract SEED file from repository
        pass
    
    def get_hmm_profile(self, family_id):
        # Extract HMM file from repository
        pass
```

### HHblits Integration (NEW REQUIREMENT)
**Current Status:** No HHblits integration
**What's Needed:**
- HHblits command wrapper
- Database/input file preparation
- Result parsing (M8 format)
- Batch job management (via SLURM)
- E-value filtering and ranking

**Dependencies:**
- HMMER tools (hmmsearch) already used in query_paperblast.py
- Similar subprocess patterns can be reused

---

## 2. REUSABLE INFRASTRUCTURE

### Pipeline Architecture Pattern
**Reference:** `pfam_curation_pipeline.py` (620 lines)

**Key Patterns to Adopt:**
```python
# Progress tracking with hidden files
PROGRESS_FILES = {
    'svn_checkout': '.svn_checkout_complete',
    'seed_extraction': '.seed_extraction_complete',
    'hmm_generation': '.hmm_generation_complete',
    'hhblits_search': '.hhblits_search_complete',
    'result_aggregation': '.result_aggregation_complete'
}

# SLURM job monitoring
def wait_for_jobs(self, timeout=3600):
    """Monitor job completion via squeue"""
    
# UniProt API integration (can reuse for accession lists)
def get_proteome_accessions(self, proteome_id):
    """Fetch from UniProt REST API"""
```

### Subprocess Execution Pattern
**Reference:** `iterate.py` (1000+ lines)

```python
def run_command(cmd, check=True, capture_output=False):
    """Robust subprocess wrapper with error handling"""
    try:
        result = subprocess.run(cmd, shell=True, check=check, 
                              capture_output=capture_output, text=True)
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}")
        raise
```

### File Organization
**Reference:** Existing curation directories

```
EXTRACTION_DIRECTORY/
├── SEED/                          # SEED alignments by family
│   ├── PF00001_SEED
│   ├── PF00002_SEED
│   └── ...
├── HMM/                           # HMM profiles by family
│   ├── PF00001.hmm
│   ├── PF00002.hmm
│   └── ...
├── HHBLITS_RESULTS/              # All-against-all results
│   ├── PF00001.vs.all.m8
│   ├── PF00002.vs.all.m8
│   └── ...
├── METRICS/                       # Quality/comparison metrics
│   ├── hhblits_summary.tsv       # Top hits, E-values, coverage
│   └── comparison_stats.json
└── .progress                      # Track completion
```

---

## 3. WORKFLOW DESIGN RECOMMENDATIONS

### Stage 1: SVN Extraction
```
Input: Family accession list (PF00001, PF00002, ...)
↓
For each family:
  - SVN checkout latest version
  - Extract SEED alignment file
  - Extract HMM profile
  - Cache in local directory
Output: Local SEED/ and HMM/ directories
```

### Stage 2: HMM Validation & Preparation
```
Input: HMM directory
↓
For each HMM:
  - Verify HMM format (hmmfile tool)
  - Extract basic stats (num_seqs, model_length)
  - Create index if needed
Output: HMM quality metrics, indexed HMMs
```

### Stage 3: All-Against-All HHblits Comparisons
```
Input: HMM directory (N families)
↓
Generate job matrix (N×N comparisons)
↓
For each pair via SLURM job:
  - hhblits -i query.hhm -d database -o output.m8
  - Filter results (E-value threshold)
  - Extract top matches
Output: HHBLITS_RESULTS/ with all comparisons
```

### Stage 4: Result Aggregation & Analysis
```
Input: All HHBLITS_RESULTS files
↓
For each family:
  - Find top N hits
  - Calculate coverage statistics
  - Identify potential merges/overlaps
  - Rank by E-value and coverage
Output: Summary report (TSV format like triage file)
```

---

## 4. FILE FORMAT RECOMMENDATIONS

### SEED Alignment (Stockholm Format)
```
# STOCKHOLM 1.0
#=GF ID     PF00001
#=GF DE     Family description
#=GF SQ     100

P12345.1/10-100       MKTLIV...FTGK
Q98765.2/45-130       MKTLIV...FTGK
//
```

### HMM Profile Format
- HMMER3/f format (.hmm files)
- Binary or ASCII interchange format
- Can be validated with `hmmstat` or `hmmscan`

### All-Against-All Results (M8 Format)
```
PF00001  PF00002  95.2  150  3  0  1  150  1  150  1.2e-45  234.5
PF00001  PF00003  87.1  120  15  1  10  120  5  120  3.4e-22  89.3
...
```

### Summary Metrics File (TSV)
```
family_query  family_target  evalue  bitscore  coverage  num_matches  best_hit
PF00001       PF00002        1.2e-45  234.5    100%      1            PF00002
PF00001       PF00003        3.4e-22  89.3     85%       5            PF00003
...
```

---

## 5. KEY TOOLS & COMMANDS

### SVN Commands
```bash
svn checkout <URL>/<family> <dest>      # Initial checkout
svn update                               # Update to latest
svn list <URL>                          # List families
svn info                                # Get repository info
```

### HMMER/HHblits Commands
```bash
hmmstat <profile>                       # Profile statistics
hmmscan -o <output> <database> <query>  # Search profile
hhblits -i <query> -d <database> -o <m8_output>
hhmake -i <alignment> -o <profile>      # Build profile
```

### Pfam Existing Tools (for reference)
```bash
pfbuild                                 # Build HMM from SEED
belvu -n 80 ALIGN -o mul                # Non-redundancy filtering
create_alignment.pl -fasta input -m     # Alignment creation
pqc-overlap-rdb family                  # Check overlaps
```

---

## 6. INTEGRATION POINTS WITH EXISTING CODE

### Can Reuse From pfam_curation_pipeline.py:
1. **Accession input handling**
   - `load_accessions_from_file()`
   - `get_proteome_accessions()`
   
2. **Progress tracking**
   - `is_step_complete()` / `mark_step_complete()`
   - Hidden file marker pattern

3. **SLURM monitoring**
   - `wait_for_jobs()` using `squeue`

4. **Error handling & logging**
   - Structured logging setup
   - Subprocess error management

### Can Reuse From iterate.py:
1. **File utilities**
   - `num_seq()` - count sequences in alignment
   - `get_alignment_hash()` - detect changes

2. **Subprocess patterns**
   - `run_command()` function

3. **Quality metrics**
   - `is_flush()` - alignment quality

### Can Reference From clan_network_viz.py:
1. **MySQL database access** (if storing results)
   ```python
   import mysql.connector
   self.connection = mysql.connector.connect(**db_config)
   ```

2. **Result visualization** (if needed)
   - Network graph creation patterns
   - Interactive HTML export

---

## 7. CONFIGURATION REQUIREMENTS

### Environment Variables / Config File Needed:
```yaml
# svn_config.yaml
svn:
  repository_url: "svn+ssh://svn.pfam.org/repos/families"
  username: "curator"
  cache_dir: "/tmp/pfam_svn_cache"
  max_retries: 3

hhblits:
  database_path: "/data/hhblits/databases"
  evalue_threshold: 1e-10
  coverage_threshold: 0.5
  num_threads: 8
  slurm_partition: "gpu"
  timeout: 3600

output:
  results_dir: "./HHBLITS_RESULTS"
  metrics_file: "comparison_metrics.tsv"
  log_dir: "./logs"

progress:
  marker_dir: "./.progress_markers"
```

### Paths to Configure:
- SVN repository URL and credentials
- HMM database location (for hhblits reference database)
- Output directories
- SLURM queue/partition
- Timeout values

---

## 8. EXPECTED OUTPUTS & DELIVERABLES

### Phase 1 Outputs (SVN + Extraction):
- [ ] Local SEED/ directory with all SEED alignments
- [ ] Local HMM/ directory with all .hmm files
- [ ] Extraction log with success/failure counts
- [ ] Cache validation report

### Phase 2 Outputs (All-Against-All):
- [ ] HHBLITS_RESULTS/ directory with all .m8 files
- [ ] SLURM job logs
- [ ] Timing/performance metrics

### Phase 3 Outputs (Analysis):
- [ ] Summary metrics file (TSV format)
- [ ] Top hits ranking per family
- [ ] Statistical summary (hits distribution, E-value ranges)
- [ ] Potential merge candidates list
- [ ] Visual comparison report (optional)

---

## 9. PERFORMANCE CONSIDERATIONS

### Optimization Strategies:
1. **Batch SVN operations** - Checkout multiple families in parallel
2. **Cache HMM files** - Avoid redundant extractions
3. **SLURM job arrays** - Submit all-against-all as job array
4. **Result streaming** - Process .m8 files on-the-fly to avoid memory issues
5. **Incremental execution** - Support resuming from failed stages

### Expected Scale:
- N families in Pfam
- N×N comparisons needed
- Each HHblits search: ~1-10 seconds depending on DB size
- Total time: Hours to days depending on infrastructure

### Memory Requirements:
- SVN cache: ~1 GB per 100 families
- HMM profiles: ~1-10 MB per profile
- Results file: ~1-100 MB depending on number of hits

---

## 10. NEXT STEPS FOR IMPLEMENTATION

### Immediate Actions:
1. [ ] Confirm SVN repository URL and authentication method
2. [ ] Identify HMM database for hhblits (reference database)
3. [ ] Set up configuration file structure
4. [ ] Create SVNManager wrapper class
5. [ ] Create HHblitsRunner wrapper class

### Development Phases:
1. **Phase 1:** SVN extraction + local caching (Week 1)
2. **Phase 2:** HHblits integration with SLURM (Week 2)
3. **Phase 3:** Result parsing and aggregation (Week 3)
4. **Phase 4:** Testing and optimization (Week 4)

### Testing Strategy:
- Start with small test set (10-20 families)
- Verify SVN extraction accuracy
- Validate HHblits result format
- Check SLURM integration
- Full pipeline run on complete dataset

---

## Files to Reference
- `/home/user/PfamCuration/SCRIPTS/pfam_curation_pipeline.py` - Pipeline template
- `/home/user/PfamCuration/SCRIPTS/iterate.py` - Tool integration patterns
- `/home/user/PfamCuration/SCRIPTS/ted_build.py` - API integration example
- `/home/user/PfamCuration/SCRIPTS/clan_network_viz.py` - Database queries

## Full Analysis Document
- `/home/user/PfamCuration/CODEBASE_ANALYSIS.md` - Complete exploration results

