
 # Napoli Project – Prioritization of Non-Coding Variants (GRCh38)

## Cohort & Objective

We analyze **Napoli sequences**:

* **5 samples** including **mothers and children**, and
* **older sequences (grep)** (from Silvia's paper).

**Objective:** Prioritize **non-coding variants** (regulatory genome) using population rarity and multi-source functional evidence, producing final **GREEN-VARAN Levels (L1–L4)** for downstream interpretation.

### Step 1 — Variant calling.
done by Silvia.

### Step 2 — Add **gnomAD** Allele Frequencies

**What it adds:** `gnomAD_AF` (global), `gnomAD_AF_nfe` (non-Finnish European).
**Biological contribution:**

* **Rarity filter.** Common variants are usually benign; rare variants (e.g., AF < 0.01, < 0.001) are more suspicious in rare disease settings.
* **Population context.** gnomAD covers large, diverse cohorts, enabling robust separation of **population polymorphisms** from potential **disease alleles**.
* **Cohort relevance.** `AF_nfe` provides a subgroup lens closer to typical European ancestry, reducing misclassification from population stratification.

We want to establish the **first pass** of prioritization by rarity, preventing downstream steps from overvaluing common alleles.

---

### Step 3 — Re-run **GREEN-VARAN** (with AF)

**What it does:** Re-evaluates variant levels using **GREEN-DB** (regulatory annotations) modulated by **rarity cutoffs** from gnomAD.
**Biological contribution:** Aligns **regulatory evidence** with **population plausibility**—rare variants in regulatory elements rise in priority.

---

### Step 4 — Add **ReMM** (Regulatory Mendelian Mutation) Score

**What it adds:** `ReMM` (0–1), a probabilistic predictor of **non-coding pathogenicity**.
**Biological contribution:**

* Integrates **functional genomics** (chromatin accessibility, histone marks), **conservation**, regulatory element context, and sequence features.
* Trained specifically to distinguish **disease-relevant non-coding variants**, complementing coding-centric tools.

So, it adds **variant-level pathogenicity** likelihood rooted in **regulatory function**, beyond mere location within a regulatory region.

---

### Step 5 — Re-run **GREEN-VARAN** (AF + ReMM)

**What it does:** Prioritization reflects **concordance of rarity** and **predicted regulatory damage**.
**Biological contribution:** A rare variant with **high ReMM** in a regulatory window becomes a **high-value candidate**.

---

### Step 6 — Add **ncER** (non-coding Essential Regulation) Percentile

**What it adds:** `ncER` percentile—the **regional intolerance** to variation in non-coding space.
**Biological contribution:**

* Identifies **genomic neighborhoods** where variation is **depleted** in human populations, implying **functional constraint**.
* Complements ReMM: whereas ReMM scores a **specific variant**, ncER scores the **surrounding region’s essentiality**.

It provides **regional constraint** context—variants in highly constrained non-coding regions are more likely to have **deleterious regulatory effects**.

---

### Step 7 — Re-run **GREEN-VARAN** (AF + ReMM + ncER)

**What it does:** Integrates **rarity**, **variant-level risk**, and **regional non-coding constraint**.
**Biological contribution:** Triangulation—when all three align (rare + high ReMM + high ncER), confidence in **regulatory impact** strengthens substantially.

---

### Step 8 — Add **FATHMM-MKL (Non-coding)** Score

**What it adds:** `FATHMM_NC`, an **independent ML predictor** of non-coding deleteriousness.
**Biological contribution:**

* Uses different features/modeling than ReMM, bringing **orthogonal evidence** (regulatory annotations, conservation, sequence context).
* Agreement between **independent predictors** (ReMM & FATHMM-MKL) increases robustness.

**Why it’s important:** Reduces **model-specific bias**—concordance across predictors is a strong signal for downstream curation.

---

### Step 9 — Final **GREEN-VARAN** Run (AF + ReMM + ncER + FATHMM)

**What it does:** Produces **L1–L4** prioritization, integrating GREEN-DB with **population rarity** and **multi-source functional evidence**.
This is a cohesive, **evidence-weighted ranking** of non-coding variants for the Napoli cohort (mothers/children + historical grep), ready for expert review, segregation analysis, and functional follow-up.

---

## Expected Key Outputs

* `embryos.chr21.gnomad.vcf.gz`
* `embryos.chr21.gnomad.remm.vcf.gz`
* `embryos.chr21.gnomad.remm.ncer.vcf.gz`
* `embryos.chr21.fullannot.vcf.gz`
* `embryos.chr21.greenvaran.final.vcf.gz` (+ `.tbi`, TSV/logs)

---

## Implementation Notes (Compact)

* **Host vs Container:**

  * **HOST:** `bcftools annotate/index/view` to add gnomAD, ReMM, ncER, FATHMM and manage headers.
  * **CONTAINER:** `greenvaran smallvars` to re-run prioritization after each annotation layer.

* **Genome Build & Contigs:** Ensure **GRCh38** across resources; fix `chr21` vs `21` mismatches before Step 2.

* **Field Names:**

  * gnomAD fields may appear as `AF`/`AF_nfe` in the resource VCF; we expose them as `gnomAD_AF`/`gnomAD_AF_nfe`.
  * FATHMM field used here: `FATHMM_NC` (adjust QC grep accordingly).

* **Indexing:** Always (re)index `*.vcf.gz` (`.tbi`) after each annotation for fast random access.

* **QC “Quick Checks”:**

  * Headers: `bcftools view -h <vcf.gz> | grep -E 'gnom|ReMM|ncER|FATHMM'`
  * Spot-check counts (e.g., how many variants gain ReMM/ncER/FATHMM).

* **Reproducibility:** Record resource versions (gnomAD build, ReMM/ncER/FATHMM release dates, GREEN-VARAN version/config).

---

## Biological Summary

* **gnomAD**: Filters by **population rarity**, focusing attention on plausible disease alleles.
* **ReMM**: Estimates **variant-level regulatory pathogenicity** in non-coding space.
* **ncER**: Highlights **intolerant regions** where variation is likely deleterious.
* **FATHMM-MKL (NC)**: Adds **independent ML evidence** for non-coding impact.

**Iterative GREEN-VARAN re-runs** integrate these layers into a final, ranked list (L1–L4) of **candidate non-coding regulatory variants** for the Napoli cohort.

---

