# Paper Outline & Lit Review Plan: STR Kinship LR Testing

---

## Key Findings from Our Simulations

### 1. Parent-Child vs Full-Sibling Tests Have Fundamentally Different Error Profiles

**Parent-child (k0=0):** The test is all-or-nothing. If any single locus shows 0 IBS sharing, the combined LR collapses to exactly 0 (because the per-locus LR = k0 = 0 and the combined LR is a product). This makes the parent-child test extremely specific — false positives are essentially zero at core_13 with correct population (2/100k for Asian, 5/100k for Cauc). But it also means a single mutation or genotyping error at one locus can destroy a true positive.

**Full siblings (k0=0.25):** The LR never zeroes out (0 shared alleles gives LR=0.25, not 0). Discrimination is gradual and depends on accumulating evidence across loci. At core_13 with correct population, FP rate is 0.15-0.19% (LR>1), dropping to near zero at autosomal_29. FN rates: ~1.7-2.1% at core_13, dropping to 0.025-0.065% at autosomal_29.

**Cross-relationship confusion:** Full siblings tested under the parent-child hypothesis show 98.3% "false positive" rate at core_13 (they look like parent-child pairs). Half-siblings tested as full siblings: ~60% exceed LR>1.

### 2. Expanding from 13 to 20 CODIS Loci Provides the Biggest Discrimination Gain

Mean log10 LR for parent-child (correct pop):
| Loci Set | Mean log10 LR | Gain from core_13 |
|----------|--------------|-------------------|
| core_13 | 4.54 | — |
| identifiler_15 | 5.23 | +0.69 |
| expanded_20 | 6.99 | +2.45 |
| supplementary (23) | 8.06 | +3.52 |
| autosomal_29 | 10.13 | +5.59 |

The 13→20 jump (+2.45 log10 units) is the single biggest improvement. Going from 20→23 adds only +1.07. The practical implication: the 2017 CODIS expansion from 13 to 20 loci was the most impactful policy change for kinship testing accuracy.

### 3. Population Misassignment Inflates LRs Toward Overconfidence

Using wrong population frequencies **inflates** LRs (does not reduce them), increasing both false positive rates and apparent confidence in true positives. This is because using the wrong frequency table makes allele sharing look "rarer" than it actually is.

- Asian<>AfAm is the worst mismatch: 4-5 log10 unit inflation
- Asian<>Cauc: 3-4 log10 unit inflation
- Hispanic mismatches tend to be less extreme

**Pooled "all" frequencies:** Slightly worse than correct population (~0.5-1 log10 unit higher LRs) but dramatically better than the worst wrong-population scenario. The key insight: **pooled frequencies are the safest default when ancestry is uncertain**, because the cost of using "all" vs correct-pop is small, but the cost of guessing wrong is enormous.

**False positive inflation example (full-sibling test, core_13):**
- Correct pop: 0.15-0.19% FP
- Wrong pop (worst case): up to 0.29% FP
- Pooled "all": ~0.19% FP (close to correct)

---

## Comparison to Prior Literature

### Ge et al. (2011) — Closest Benchmark
- 13 CODIS loci only, Caucasian population only
- Included mutation model (we don't)
- Their parent-child FN at KI>=1000: 14.1% vs our ~2.5%
- Difference explained by their mutation model (mutations create apparent non-sharing at loci, killing the combined LR for parent-child)
- They did not test multiple populations or loci panels

### What's Novel in Our Work
1. **Systematic 5-panel comparison** (core_13, identifiler_15, expanded_20, supplementary, autosomal_29) — no prior study has done this
2. **Full factorial population mismatch analysis** (4 populations x 5 tested-population assumptions x 5 loci sets x 6 relationships x 2 tested hypotheses) — unprecedented scale
3. **Cross-relationship confusion matrix** — systematic quantification of how often one relationship is confused for another
4. **Pooled vs population-specific frequencies** — direct empirical comparison showing pooled is safer than wrong-specific

---

## Lit Review Needs (for RA)

### Category 1: Kinship LR Methodology & Validation
**Goal:** Establish the standard framework and identify what simulation studies have been done.
- Search: "likelihood ratio kinship testing STR validation simulation"
- Search: "kinship index STR forensic genetics false positive"
- Key question: Who else has done large-scale simulation studies? At what scale? Which loci panels?
- Look for: Ge et al. (2011), Kling et al., Buckleton et al., Gill et al.

### Category 2: CODIS Expansion Impact
**Goal:** Frame the 2017 expansion from 13->20 core loci and what's known about its effect on kinship testing.
- Search: "CODIS expansion 20 loci kinship discrimination power"
- Search: "expanded CODIS STR panel forensic identification"
- Key question: Has anyone quantified the kinship testing improvement from the expansion? (We think not systematically)

### Category 3: Population Allele Frequency Effects on LR
**Goal:** Document what's known about population misassignment and LR bias.
- Search: "population allele frequency mismatch likelihood ratio forensic"
- Search: "reference population selection STR kinship bias"
- Search: "pooled allele frequencies forensic kinship testing"
- Key question: Has anyone shown that wrong-population frequencies inflate LRs? Has anyone recommended pooled frequencies as a safer default?

### Category 4: Error Profiles of Different Kinship Tests
**Goal:** Understand the theoretical basis for why parent-child and full-sibling tests behave so differently.
- Search: "parent child sibling kinship test error rate comparison STR"
- Search: "IBS sharing probability unrelated individuals STR"
- Key question: Is the all-or-nothing behavior of parent-child testing (k0=0) well-documented in the literature? Is it discussed as a practical concern?

### Category 5: Mutation Effects on Kinship LR
**Goal:** Contextualize our no-mutation model and explain the difference from Ge et al.
- Search: "STR mutation rate kinship likelihood ratio effect"
- Search: "mutation model parent offspring testing false negative"
- Key question: How do mutations affect parent-child vs full-sibling tests differently? This is important because our FN rates are lower than Ge et al.'s, and we need to explain why.

---

## Proposed Paper Structure

### Title (draft)
"Systematic evaluation of STR-based kinship likelihood ratios across loci panels, populations, and reference frequency assumptions"

### Abstract
[Summary of 3 main findings + scale of simulation]

### Introduction
- Forensic kinship testing using STR likelihood ratios is standard practice
- The field has moved from 13 to 20+ CODIS loci, but systematic evaluation of the impact on kinship testing is limited
- Choice of reference population frequencies is known to matter, but the direction and magnitude of bias from misassignment is not well characterized
- We present a large-scale simulation study evaluating...

### Methods
- Simulation framework (modules 1-6): genotype simulation under HWE, IBD-based relative generation, LR calculation
- 6 relationships x 4 populations x 5 loci panels x 5 reference-population assumptions
- 100k unrelated pairs + 20k related pairs per combination
- No mutation model (discuss as limitation, compare to Ge et al.)
- FPR-based cutoff calculation methodology

### Results
1. **Error profiles differ by relationship type** (parent-child all-or-nothing vs full-sibling gradual)
2. **Loci panel comparison** (13->20 biggest jump, diminishing returns after)
3. **Population mismatch inflates LRs** (direction, magnitude, worst-case pairs)
4. **Pooled frequencies as safest default** (small cost vs correct, large benefit vs wrong)

### Discussion
- Practical implications for forensic labs choosing panels and reference populations
- The case for pooled frequencies when ancestry is uncertain
- Limitations: no mutation model, HWE assumption, no population substructure
- Comparison to Ge et al. and other prior work

### Figures (planned)
1. FP rate heatmap by loci set x population (PC vs FS)
2. TP rates across thresholds and loci sets
3. Cross-relationship confusion matrix
4. Population mismatch effect on LR magnitude
5. Correct vs pooled vs wrong-population FP comparison

---

## Pending Analyses

### Focal Simulation (modules 7-8)
Modules exist but are not yet in the full pipeline. The focal approach is more realistic:
- Generate focal individual with known family structure
- Generate pool of unrelated individuals
- Test focal against each pool member AND known relatives
- Better mimics real forensic database searches

This would strengthen the paper by adding a "database search" scenario alongside the current pairwise analysis.
