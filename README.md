
 # Deliver - Prioritization of Non-Coding Variants

This repository contains scripts and results related to the non-coding prioritization of the DELIVER project from the CNR (PRIN prot. 2020J84FAM).
## Cohort & Objective


We analyze **Napoli sequences**:

* **5 samples** including **mothers and children**, and
* **older sequences (grep)** from "Prioritization of putatively detrimental variants in euploid miscarriages" (https://www.nature.com/articles/s41598-022-05737-3).

**Objective:** Prioritize **non-coding variants** (regulatory genome) using population rarity and multi-source functional evidence, producing final **GREEN-VARAN Levels (L1–L4)** for downstream interpretation.

### Step 1 — Variant calling.
Already done. 

### Step 2 — Add **gnomAD** Allele Frequencies

Source: https://gnomad.broadinstitute.org/

**What it adds:** `gnomAD_AF` (global), `gnomAD_AF_nfe` (non-Finnish European).
**Biological contribution:**

* **Rarity filter.** Common variants are usually benign; rare variants (e.g., AF < 0.01, < 0.001) are more suspicious in rare disease settings.
* **Cohort** `AF_nfe` subgroup for this population (not for the Malaysia case).

We want to establish the **first pass** of prioritization by rarity, preventing downstream steps from overvaluing too many alleles.

---

### Step 3 — Run **GREEN-VARAN** (with AF)

Source: https://pubmed.ncbi.nlm.nih.gov/35234913/

**What it does:** Re-evaluates variant levels using **GREEN-DB** (regulatory annotations) modulated by **rarity cutoffs** from gnomAD.
**Biology:** Aligns **regulatory evidence** with **population plausibility**—rare variants in regulatory elements rise in priority.

---

### Step 4 — Add **ReMM** (Regulatory Mendelian Mutation) Score

Source: https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giad024/7135629?login=false

**What it adds:** `ReMM` (0–1), a probabilistic predictor of **non-coding pathogenicity**.
**Biological contribution:**

* Integrates **functional genomics** (chromatin accessibility, histone marks), **conservation**, regulatory element context, and sequence features.
* Trained specifically to distinguish **disease-relevant non-coding variants**, complementing coding-centric tools.

So, it adds **variant-level pathogenicity** likelihood rooted in **regulatory function**, beyond mere location within a regulatory region.

---

### Step 5 — Re-run **GREEN-VARAN** (AF + ReMM)

**What it does:** Prioritization reflects **concordance of rarity** and **predicted regulatory damage**.
**Biology:** A rare variant with **high ReMM** in a regulatory window becomes a **high-value candidate**.

---

### Step 6 — Add **ncER** (non-coding Essential Regulation) Percentile

Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC6868241/

**What it adds:** `ncER` percentile—the **regional intolerance** to variation in non-coding space.
**Biology:**

* Identifies **genomic neighborhoods** where variation is **depleted** in human populations, implying **functional constraint**.
* Complements ReMM: whereas ReMM scores a **specific variant**, ncER scores the **surrounding region’s essentiality**.

It provides **regional constraint** context—variants in highly constrained non-coding regions are more likely to have **deleterious regulatory effects**.

---

### Step 7 — Run **GREEN-VARAN** again (AF + ReMM + ncER)

**What it does:** Integrates **rarity**, **variant-level risk**, and **regional non-coding constraint**.
**Biological contribution:** Triangulation, when all three align (rare + high ReMM + high ncER), confidence in **regulatory impact** strengthens.

---

### Step 8 — Add **FATHMM-MKL (Non-coding)** Score

Source: https://academic.oup.com/bioinformatics/article/34/3/511/4104409

**What it adds:** `FATHMM_NC`, an **independent ML predictor** of non-coding deleteriousness.
**Biological contribution:**

* Uses different features/modeling than ReMM, bringing **orthogonal evidence** (regulatory annotations, conservation, sequence context).
* Agreement between **independent predictors** (ReMM & FATHMM-MKL) increases robustness.

**Why?:** Reduces **model-specific bias**, concordance across predictors is a strong signal for downstream curation.

---

### Step 9 — Final **GREEN-VARAN** Run (AF + ReMM + ncER + FATHMM)

It produces **L1–L4** prioritization, integrating GREEN-DB with **population rarity** and **multi-source functional evidence**.
This is a cohesive, **evidence-weighted ranking** of non-coding variants for the Napoli cohort (mothers/children + historical grep), ready for expert review, segregation analysis, and functional follow-up.

---

Finally, de novo mutations can be identified by comparing with mothers data.

## Expected outputs

* `embryos.chr.gnomad.vcf.gz`
* `embryos.chr.gnomad.remm.vcf.gz`
* `embryos.chr.gnomad.remm.ncer.vcf.gz`
* `embryos.chr.fullannot.vcf.gz`
* `embryos.chr.greenvaran.final.vcf.gz` (+ `.tbi`, TSV/logs)

---

