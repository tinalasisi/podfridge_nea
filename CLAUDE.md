# CLAUDE.md — Project Context for Claude Code Sessions

## Project Overview

**podfridge_nea** is a forensic genetics simulation pipeline that:
1. Simulates STR (Short Tandem Repeat) genotype pairs at various kinship relationships
2. Computes likelihood ratios (LRs) for kinship hypotheses
3. Evaluates false positive/negative rates across populations, loci panels, and population assumptions

The goal is a forensic genetics paper evaluating STR-based kinship LR testing.

## Repository Structure

### Data Files (`data/`)
- `df_allelefreq_combined.csv` — Allele frequency table (columns: population, marker, allele, frequency). Populations: AfAm, Asian, Cauc, Hispanic, all (pooled)
- `core_CODIS_loci.csv` — Loci set membership (columns: locus, core_13, identifiler_15, expanded_20, supplementary)
- `kinship_coefficients.csv` — IBD coefficients (k0, k1, k2) for 6 relationships: parent_child (0,1,0), full_siblings (0.25,0.5,0.25), half_siblings (0.5,0.5,0), cousins (0.875,0.125,0), second_cousins (0.9375,0.0625,0), unrelated (1,0,0)
- `focal_family_configs.csv` — Family structure configs for module 7

### Code Modules (`code/`)

**Core simulation pipeline (modules 1-6):**
- `LR_kinship_utility_functions.R` — Shared constants (FALLBACK_FREQ = 5/(2*1036)), allele frequency loading, loci set definitions, LR formula, genotype matching. This is sourced by almost everything.
- `module1_allele_simulator.R` — Simulate single allele from frequency distribution
- `module2_STR_profile_simulator.R` — Simulate full STR profile (all loci) under HWE
- `module3_related_individual_simulator.R` — Simulate related individual given focal profile + IBD coefficients
- `module4_single_locus_LR.R` — Calculate single-locus LR for a pair
- `module5_combined_LR.R` — Multiply single-locus LRs across loci sets to get combined LR
- `module6_single_combo_pair_generator.R` — Wrapper: generate batch of pairs for one population x relationship combo

**Focal simulation (modules 7-8, NOT YET USED IN FULL PIPELINE):**
- `module7_single_pop_focal_generator.R` — Generate focal individuals with flexible family structures
- `module8_unrelated_pool_generator.R` — Generate pool of unrelated individuals

**Analysis and statistics:**
- `module9_combinedLR_stats_functions.R` — Summary stats, FPR-based cutoff calculation, proportion exceeding cutoffs
- `analyze_lr_outputs.R` — Main post-processing: reads combined_LR files, computes stats per population
- `sanity_check_analysis.R` — Comprehensive FP/FN rates across ALL combinations
- `sanity_check_figures.R` — 6 diagnostic figures from sanity check output

**SLURM submission scripts:**
- `sim_pairs.sh` — Submit pair simulation (20 populations x 5 relationships x 10 chunks of 2000 = 1000 jobs)
- `sim_pairs_unrelated.sh` — Submit unrelated pair simulation (4 pops x 10 chunks of 10000 = 40 jobs)
- `lr_submission.sh` / `lr_wrapper.sh` / `lr_wrapper.R` — Single-locus LR calculation
- `combined_lr_submission.sh` / `combined_lr.sh` / `combined_lr_wrapper.R` — Combined LR calculation
- `analyze_lr_outputs.sh` — Analysis submission
- `sanity_check.sh` — Sanity check submission
- `run_pipeline_bare.sh` — Step-by-step pipeline commands

**Plotting:**
- `plots_matched.R` — Correct-population LR distributions
- `plots_mismatched.R` — Population mismatch LR distributions
- `plots_proportion_exceeding_cutoffs.R` — Proportion exceeding FPR-based thresholds

### Pipeline Flow
```
sim_pairs.R → output/pairs_*.csv
  → lr_wrapper.R → output/LR/LR_*.csv
    → combined_lr_wrapper.R → output/combined_LR/combined_LR_*.csv
      → analyze_lr_outputs.R → output/lr_analysis_YYYYMMDD/
      → sanity_check_analysis.R → output/sanity_check/
```

## Key Technical Details

### LR Formula
The core LR calculation in `LR_kinship_utility_functions.R`:
- `shared_alleles == 0`: returns `k0` (for parent_child where k0=0, this zeros out the combined LR)
- `shared_alleles == 1`: returns `k0 + k1/Rxp` where Rxp depends on genotype match
- `shared_alleles == 2`: returns `k0 + k1/Rxp + k2/Rxu`

**Critical insight**: For parent-child (k0=0), any single locus with 0 IBS sharing makes the entire combined LR = 0 (product of per-locus LRs). For full siblings (k0=0.25), LR never zeroes out because 0 shared alleles still gives LR=0.25.

### FALLBACK_FREQ
`5 / (2 * 1036) ≈ 0.002414` — Used when an allele has zero frequency in a population's table.

### Simulation Scale (current run)
- 100,000 unrelated pairs per population (4 populations)
- 20,000 related pairs per population x relationship (4 populations x 5 relationships)
- 5 loci sets: core_13, identifiler_15, expanded_20, supplementary (23 loci), autosomal_29
- 5 population frequency tables tested: AfAm, Asian, Cauc, Hispanic, all

### FPR-Based Cutoffs
Computed from quantiles of the unrelated pairs' LR distribution at 1%, 0.1%, 0.01% FPR. For full siblings, the 1% FPR cutoff is ~0.09 (below LR=1), because 99% of unrelated pairs already have LR < 0.09 for the full-sibling hypothesis. This is correct behavior, not a bug.

## Completed Work & Key Findings

### Sanity Check Results (January 2026 run)
All patterns match population genetics theory:
1. **Parent-child FP**: Essentially zero at core_13 with correct population (2/100k for Asian, 5/100k for Cauc)
2. **Full-sibling FP**: 0.15-0.19% at core_13 correct pop (LR>1), drops rapidly with more loci
3. **Population mismatch inflates LRs**: Using wrong population frequencies INFLATES LRs, increasing both FP rates and apparent TP confidence. Asian↔AfAm is the worst mismatch (4-5 log10 units inflation)
4. **Pooled "all" frequencies**: Slightly worse than correct population but dramatically better than worst wrong population. Safest bet when population is unknown.
5. **Expanding from 13→20 loci provides the biggest improvement** in discrimination power

### Paper Story (3 main points)
1. Parent-child vs full-sibling tests have fundamentally asymmetric error profiles (all-or-nothing vs gradual)
2. Expanding from 13 to 20 CODIS loci provides the biggest gain; going beyond 20 has diminishing returns
3. Population misassignment inflates LRs toward overconfidence — pooled frequencies are the safest choice when ancestry is uncertain

### Literature Context
- **Ge et al. (2011)** is the closest prior work: 13 CODIS loci, Caucasian only, with mutations. Their parent-child FN at KI≥1000 was 14.1% vs our ~2.5% (difference due to their mutation model)
- **Novel contributions**: Systematic 5-panel comparison, full factorial population mismatch analysis, cross-relationship confusion matrix at scale

## Pending Work

### Focal Simulation (Priority)
Modules 7 and 8 exist but are **not yet integrated into the full pipeline**. The focal simulation approach:
- Generate a "focal" individual with known relatives (using module 7)
- Generate a pool of unrelated individuals (using module 8)
- Test the focal against each pool member AND their known relatives
- This is more realistic than the current pairwise approach (which tests pre-defined pairs)

### Other Pending
- Lit review framing for RA (search queries were drafted but not finalized)
- Formal paper writing

## SLURM Environment
- Account: `tlasisi0`, partition: `standard`
- R module: `module load Rtidyverse`
- Email notifications: `nicolead@umich.edu`
- Typical resources: 96G RAM, 2hr time limit for analysis; 4G RAM per simulation job

## Git Conventions
- Branch naming: `claude/<random-words>` for Claude Code sessions
- PRs go to `main` on `tinalasisi/podfridge_nea` (fork of `NicoleAdams-sci/podfridge_nea`)
