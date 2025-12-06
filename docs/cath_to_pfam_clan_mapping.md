# CATH to Pfam Clan Mapping

This document describes the process for creating a mapping between CATH structural classifications and Pfam clans. The mapping is based on domain overlaps on UniProt sequences.

## Purpose

When curating new Pfam families from TED (The Encyclopedia of Domains) data, TED domains often have CATH classifications assigned (e.g., `3.40.50.300`). This mapping allows us to automatically suggest appropriate Pfam clans for new families based on their CATH classification.

## Data Sources

| Data | Source | Description |
|------|--------|-------------|
| TED domain summary | [Zenodo](https://zenodo.org/records/13908086) | `ted_365m.domain_summary.cath.globularity.taxid.tsv.gz` (~20GB) |
| Pfam domain regions | `pfam_live` database | Table `pfamA_reg_full_significant` |
| Pfam clan membership | `pfam_live` database | Tables `clan_membership` and `clan` |

## Prerequisites

- Access to `pfam_live` MySQL database with `~/.my.cnf` configured
- TED domain summary file downloaded from Zenodo
- ~50GB disk space for intermediate files
- Python 3 with standard library

## Step-by-Step Process

All commands assume you are working in `/nfs/production/agb/pfam/data/TED/`

### Step 1: Download TED Domain Summary (if not already present)

Download from Zenodo: https://zenodo.org/records/13908086

File: `ted_365m.domain_summary.cath.globularity.taxid.tsv.gz` (~20GB)

### Step 2: Extract TED Domains with CATH Labels

Extract only high-confidence domains that have a CATH classification. The chopping column is preserved to enable segment-aware overlap calculation for multi-segment domains:

```bash
zcat ted_365m.domain_summary.cath.globularity.taxid.tsv.gz | awk -F'\t' '
$3=="high" && $14!="-" {
    split($1, a, "-"); uniprot = a[2]
    n = split($1, b, "_"); ted_suffix = b[n]
    print uniprot "\t" ted_suffix "\t" $4 "\t" $14 "\t" $15
}' | gzip > ted_cath_high_segments.tsv.gz
```

**Input columns from TED file:**
- Column 1: TED ID (e.g., `AF-A0A000-F1-model_v4_TED01`)
- Column 3: Consensus level (`high`, `medium`)
- Column 4: Chopping/boundaries (e.g., `11-41_290-389` for multi-segment, `100-250` for single segment)
- Column 14: CATH label (e.g., `3.40.50.300` or `-`)
- Column 15: CATH assignment level (`H` for homologous superfamily, `T` for topology/fold)

**Output format:** `uniprot_acc, ted_suffix, chopping, cath_label, cath_level`

**Expected output:** ~174 million rows, ~1.3GB compressed

### Step 3: Extract Pfam Clan Membership

```bash
mysql --defaults-file=~/.my.cnf pfam_live --quick --skip-column-names -e "
SELECT m.pfamA_acc, m.clan_acc, c.clan_id
FROM clan_membership m
JOIN clan c ON m.clan_acc = c.clan_acc
" > pfam_clan_membership.tsv
```

**Output format:** `pfamA_acc, clan_acc, clan_id`

**Expected output:** ~13,600 rows

### Step 4: Extract Pfam Domain Regions for Clan Families

Only extract regions for families that belong to a clan:

```bash
mysql --defaults-file=~/.my.cnf pfam_live --quick --skip-column-names -e "
SELECT r.pfamseq_acc, r.seq_start, r.seq_end, r.pfamA_acc
FROM pfamA_reg_full_significant r
WHERE r.in_full = 1
AND r.pfamA_acc IN (SELECT pfamA_acc FROM clan_membership)
" | gzip > pfam_clan_regions.tsv.gz
```

**Output format:** `pfamseq_acc, seq_start, seq_end, pfamA_acc`

**Expected output:** ~144 million rows

### Step 5: Sort Both Files by UniProt Accession

The mapping script uses a merge-join algorithm that requires both files to be sorted:

```bash
# Sort TED file
zcat ted_cath_high_segments.tsv.gz | sort -S 1G -k1,1 | gzip > ted_cath_high_segments_sorted.tsv.gz

# Sort Pfam file
zcat pfam_clan_regions.tsv.gz | sort -S 1G -k1,1 | gzip > pfam_clan_regions_sorted.tsv.gz
```

**Note:** Sorting may take 20-30 minutes per file. The `-S 1G` flag limits memory usage.

### Step 6: Run the Mapping Script

```bash
python3 /path/to/PfamCuration/SCRIPTS/cath_to_pfam_mapping.py \
    --ted-file ted_cath_high_segments_sorted.tsv.gz \
    --pfam-file pfam_clan_regions_sorted.tsv.gz \
    --clan-file pfam_clan_membership.tsv \
    --output-dir .
```

**Options:**
- `--min-overlap`: Minimum overlap fraction (default: 0.5 = 50% in both directions)
- `--progress-interval`: Print progress every N accessions (default: 1,000,000)

## Output Files

### cath_to_pfam_clan_mapping.tsv

Main mapping file with confidence scores:

| Column | Description |
|--------|-------------|
| `cath_label` | CATH classification (3 or 4 levels) |
| `cath_level` | `H` (homologous superfamily) or `T` (topology/fold) |
| `clan_acc` | Pfam clan accession (e.g., `CL0063`) |
| `clan_id` | Pfam clan name (e.g., `FAD_DHS`) |
| `overlap_count` | Number of domain pair overlaps supporting this mapping |
| `confidence` | Fraction of this CATH label's overlaps that map to this clan |
| `top_families` | Top 3 Pfam families contributing to this mapping |

**Example:**
```tsv
cath_label      cath_level  clan_acc  clan_id    overlap_count  confidence  top_families
3.40.50.300     H           CL0063    FAD_DHS    12456          0.92        PF00070,PF00175,PF00890
1.10.10         T           CL0015    Homeobox   5621           0.75        PF00046,PF00520
```

### cath_to_pfam_family_counts.tsv

Detailed per-family counts for analysis:

| Column | Description |
|--------|-------------|
| `cath_label` | CATH classification |
| `cath_level` | `H` or `T` |
| `pfamA_acc` | Pfam family accession |
| `count` | Number of overlapping domain pairs |

## Algorithm Details

### Domain Overlap Calculation

Two domains are considered overlapping if **both** of the following are true:
- ≥50% of the TED domain is covered by the Pfam domain
- ≥50% of the Pfam domain is covered by the TED domain

This bidirectional requirement ensures meaningful structural correspondence.

### Multi-Segment Domain Handling

TED domains can be discontinuous (e.g., `11-41_290-389` for a domain with two segments). For these multi-segment domains:

1. The chopping string is preserved in the data file
2. Overlap is calculated **per segment** and summed
3. The total overlap is compared against the total TED domain length (sum of segment lengths)

This prevents false positives from nested domain scenarios where a contiguous interpretation would incorrectly inflate overlaps.

### Merge-Join Streaming

The script processes both sorted files simultaneously without loading them entirely into memory:

1. Read lines from both files
2. Compare UniProt accessions
3. Advance the file with the "smaller" accession
4. When accessions match, collect all domains for that protein and find overlaps
5. Accumulate counts and move to next accession

This allows processing hundreds of millions of domains with minimal memory (~few GB for the counts dictionary).

### CATH Classification Levels

- **H-level** (4 numbers, e.g., `3.40.50.300`): Homologous superfamily - proteins share a common ancestor
- **T-level** (3 numbers, e.g., `3.40.50`): Topology/fold - proteins share the same fold but may not be homologous

Both levels are included in the mapping. H-level mappings are more precise; T-level provides broader coverage.

## Updating the Mapping

The mapping should be regenerated when:
- New Pfam release adds/modifies clan membership
- New TED release with updated CATH assignments

Run Steps 3-6 to regenerate with updated data.

## File Sizes (Reference)

| File | Rows | Compressed Size |
|------|------|-----------------|
| `ted_365m.domain_summary.cath.globularity.taxid.tsv.gz` | 365M | ~20GB |
| `ted_cath_high_segments.tsv.gz` | 174M | ~1.3GB |
| `pfam_clan_regions.tsv.gz` | 144M | ~2GB |
| `pfam_clan_membership.tsv` | 13.6K | ~500KB |

## Contact

Created: December 2024
Author: Alex Bateman
