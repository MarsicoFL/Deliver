#!/usr/bin/env bash
set -euo pipefail

# ===== Parameters =====
CHR=20
SIF="/lustre/home/enza/franco/green-varan/greenvaran_1.3.3.sif"

# ===== Resources (HOST) =====
RES="/lustre/home/enza/franco/green-varan/resources"
RES38="$RES/GRCh38"

GNOMAD_VCF="$RES38/GRCh38_gnomad.genomes.chr${CHR}.vcf.gz"   # per-chromosome
GNOMAD_TBI="${GNOMAD_VCF}.tbi"

REMM_TSV="$RES38/GRCh38_ReMM.tsv.gz"
NCER_BED="$RES38/GRCh38_ncER_perc.bed.gz"
GREEN_DB_HOST="$RES38/GRCh38_GREEN-DB.bed.gz"

# ===== Paths INSIDE the container =====
GREEN_DB_CNTR="/opt/green-varan/resources/GRCh38/GRCh38_GREEN-DB.bed.gz"
DB_SCHEMA="/opt/green-varan/config/greendb_schema_v2.5.json"
GV_CONF="/opt/green-varan/config/prioritize_smallvars.json"

# ===== IO =====
IN_DIR="/lustre/home/enza/franco/deliver/samples/grep/"
OUT_DIR="/lustre/home/enza/franco/deliver/noncod_prior"
mkdir -p "$OUT_DIR"

IN_VCF="$IN_DIR/grepEmbryos.chr${CHR}.gatk.vcf.gz"

VCF_GNOMAD="$OUT_DIR/embryos.chr${CHR}.gnomad.vcf.gz"
VCF_GNOMAD_REMM="$OUT_DIR/embryos.chr${CHR}.gnomad.remm.vcf.gz"
VCF_GNOMAD_REMM_NCER="$OUT_DIR/embryos.chr${CHR}.gnomad.remm.ncer.vcf.gz"
VCF_FULL="$OUT_DIR/embryos.chr${CHR}.fullannot.vcf.gz"

VCF_GV_STEP="$OUT_DIR/embryos.chr${CHR}.greenvaran.vcf.gz"
VCF_GV_FINAL="$OUT_DIR/embryos.chr${CHR}.greenvaran.final.vcf.gz"

# ===== Pre-checks =====
[[ -s "$GNOMAD_VCF" ]] || { echo "ERROR: Missing $GNOMAD_VCF"; exit 2; }
[[ -s "$GNOMAD_TBI" ]] || { echo "ERROR: Missing index $GNOMAD_TBI"; exit 2; }
[[ -s "$IN_VCF" ]]     || { echo "ERROR: Missing input $IN_VCF"; exit 2; }
[[ -s "$GREEN_DB_HOST" ]] || { echo "ERROR: Missing $GREEN_DB_HOST"; exit 2; }


# ===== STEP 2 — Add gnomAD AF (HOST) =====
echo "[STEP 2] Adding gnomAD AF from $(basename "$GNOMAD_VCF")…"
bcftools annotate \
  -a "$GNOMAD_VCF" \
  -c CHROM,POS,REF,ALT,INFO/gnomAD_AF:=INFO/AF,INFO/gnomAD_AF_nfe:=INFO/AF_nfe \
  "$IN_VCF" \
  -Oz -o "$VCF_GNOMAD"
bcftools index -f "$VCF_GNOMAD"
bcftools view -h "$VCF_GNOMAD" | grep -E 'gnomAD_AF|gnomAD_AF_nfe' || true

# ===== STEP 3 — GREEN-VARAN with gnomAD (CONTAINER) =====
# NOTE: Use container path for --db to avoid host-path open errors.
echo "[STEP 3] GREEN-VARAN with gnomAD…"
apptainer exec -B "$RES:/opt/green-varan/resources" -B "/lustre/home/enza/franco/deliver:/deliver" "$SIF" \
  greenvaran smallvars \
  --invcf "/deliver/noncod_prior/embryos.chr${CHR}.gnomad.vcf.gz" \
  --outvcf "/deliver/noncod_prior/embryos.chr${CHR}.greenvaran.vcf.gz" \
  --db "$GREEN_DB_CNTR" --dbschema "$DB_SCHEMA" --config "$GV_CONF"

# ===== STEP 4 — Add ReMM (HOST) =====
echo "[STEP 4] Adding ReMM…"
bcftools annotate \
  -a "$REMM_TSV" \
  -h <(echo '##INFO=<ID=ReMM,Number=1,Type=Float,Description="ReMM score (GRCh38)">') \
  -c CHROM,POS,ReMM \
  "$VCF_GNOMAD" \
  -Oz -o "$VCF_GNOMAD_REMM"
bcftools index -f "$VCF_GNOMAD_REMM"

# ===== STEP 5 — GREEN-VARAN with ReMM (CONTAINER) =====
# Ensure --db uses the container path again.
echo "[STEP 5] GREEN-VARAN with ReMM…"
apptainer exec -B "$RES:/opt/green-varan/resources" -B "/lustre/home/enza/franco/deliver:/deliver" "$SIF" \
  greenvaran smallvars \
  --invcf "/deliver/noncod_prior/embryos.chr${CHR}.gnomad.remm.vcf.gz" \
  --outvcf "/deliver/noncod_prior/embryos.chr${CHR}.greenvaran.vcf.gz" \
  --db "$GREEN_DB_CNTR" --dbschema "$DB_SCHEMA" --config "$GV_CONF"

# ===== STEP 6 — Add ncER (HOST) =====
echo "[STEP 6] Adding ncER…"
HDR_NCER="$OUT_DIR/hdr.ncER.txt"
echo '##INFO=<ID=ncER,Number=1,Type=Float,Description="ncER percentile (noncoding intolerance)">' > "$HDR_NCER"

bcftools annotate \
  -h "$HDR_NCER" -a "$NCER_BED" -c CHROM,FROM,TO,INFO/ncER \
  "$VCF_GNOMAD_REMM" \
  -Oz -o "$VCF_GNOMAD_REMM_NCER"
bcftools index -f "$VCF_GNOMAD_REMM_NCER"

# ===== STEP 7 — GREEN-VARAN with ncER (CONTAINER) =====
# Again, use container path for --db.
echo "[STEP 7] GREEN-VARAN with ncER…"
apptainer exec -B "$RES:/opt/green-varan/resources" -B "/lustre/home/enza/franco/deliver:/deliver" "$SIF" \
  greenvaran smallvars \
  --invcf "/deliver/noncod_prior/embryos.chr${CHR}.gnomad.remm.ncer.vcf.gz" \
  --outvcf "/deliver/noncod_prior/embryos.chr${CHR}.greenvaran.vcf.gz" \
  --db "$GREEN_DB_CNTR" --dbschema "$DB_SCHEMA" --config "$GV_CONF"

# ===== STEP 8 — Add FATHMM-MKL_NC (HOST) =====
echo "[STEP 8] Adding FATHMM-MKL_NC…"
bcftools annotate \
  -a "$RES38/GRCh38_FATHMM-MKL_NC.tsv.gz" \
  -h <(echo '##INFO=<ID=FATHMM_NC,Number=1,Type=Float,Description="FATHMM-MKL non-coding score (GRCh38)">') \
  -c CHROM,POS,REF,ALT,FATHMM_NC \
  "$VCF_GNOMAD_REMM_NCER" \
  -Oz -o "$VCF_FULL"
bcftools index -f "$VCF_FULL"

# ===== STEP 9 — GREEN-VARAN final (CONTAINER) =====
# Final call also fixed to use container path for --db.
echo "[STEP 9] GREEN-VARAN FINAL…"
apptainer exec -B "$RES:/opt/green-varan/resources" -B "/lustre/home/enza/franco/deliver:/deliver" "$SIF" \
  greenvaran smallvars \
  --invcf "/deliver/noncod_prior/embryos.chr${CHR}.fullannot.vcf.gz" \
  --outvcf "/deliver/noncod_prior/embryos.chr${CHR}.greenvaran.final.vcf.gz" \
  --db "$GREEN_DB_CNTR" --dbschema "$DB_SCHEMA" --config "$GV_CONF"

echo "Done ✅  Final output: $VCF_GV_FINAL"
