# Pfam Curator Skill v1.1

## Overview

This skill guides Claude in creating and editing protein family entries for the Pfam database, following proper curation standards.

## Key Features

- **Correct DESC file ordering**: References BEFORE CC lines (this was the critical fix)
- **Factual annotation standards**: Emphasis on experimental evidence only
- **British English guidelines**: Proper spellings throughout
- **Gene names over DUFs**: Guidance to use actual gene names when available
- **Domain vs protein function distinction**: Clear guidelines on what to attribute to domains

## Files Included

### SKILL.md
Main skill file with:
- Core principles
- Field order (with complete CBS domain example)
- Naming conventions
- Quick reference commands
- CC line writing guidelines
- Type field selection
- Clan assignment methods

### references/desc-format.md
Complete DESC field specifications with:
- **CRITICAL field order section at top with example**
- Detailed format for each field type
- Header, Reference, and Comment section specifications
- Automated vs manual field distinctions

### references/writing-guidelines.md
Annotation writing standards:
- British English spelling guide
- What to include/avoid in annotations
- Domain vs protein function guidelines
- Citation standards
- Quality checklist

### references/workflows.md
Step-by-step procedures for:
- Building new families
- Iterating existing families
- Editing annotations
- Adding to clans
- Handling overlaps
- Quick reference commands

## Critical Fix in v1.1

**Problem**: Previous version showed incorrect field order with CC lines before references

**Solution**: 
- Added explicit "CRITICAL: Field Order" sections in both SKILL.md and desc-format.md
- Included complete example (CBS domain) showing correct order
- Added multiple warnings throughout documentation
- Shows three distinct sections: Header → Reference → Comment

**Correct order**:
1. Header section (ID through TP, optional PI/CL)
2. Reference section (WK, RC, RN, RM, RT, RA, RL, DR) - **BEFORE CC**
3. Comment section (CC lines) - **ALWAYS LAST**

## Installation

For Claude Desktop App:
1. Download the pfam-curator.zip file
2. In Claude Desktop, go to Settings > Skills
3. Click "Add Skill" and select the downloaded file

For Claude.ai:
- Upload the skill through the Projects interface

## Usage

The skill activates when:
- Creating or editing Pfam family DESC files
- Annotating protein domains and families
- Assigning families to clans
- Following Pfam ID and naming conventions
- Writing CC line annotations

## Support

For questions about Pfam curation, refer to:
- Official Pfam documentation
- Your team's experienced curators
- The skill's reference files for detailed guidance
