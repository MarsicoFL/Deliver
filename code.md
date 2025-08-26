CONTEXT
- Sample VCF (chr21):  /lustre/home/enza/franco/deliver/samples/napoli/embryos/embryos.chr21.gatk.vcf.gz
- Work/output dir:      /lustre/home/enza/franco/deliver/noncod_prior
- GREEN-VARAN image:    /lustre/home/enza/franco/green-varan/greenvaran_1.3.3.sif
- Resources (GRCh38):   /lustre/home/enza/franco/green-varan/resources/GRCh38/
- You already ran Level-0 and have: embryos.chr21.greenvaran.vcf.gz (no AF/scores yet)

ENV RULE
- [HOST] = run on HPC shell (your local bcftools).
- [INSIDE CONTAINER] = run via Apptainer; only greenvaran here.

──────────────────────────────────────────────────────────────────────────────

STEP 2 — Add gnomAD AF to the sample VCF  [HOST, one-liner]
Adds: INFO/gnomAD_AF, INFO/gnomAD_AF_nfe

bcftools annotate \
  -a /lustre/home/enza/franco/green-varan/resources/GRCh38/GRCh38_gnomad.genomes.fixed.vcf.gz \
  -c CHROM,POS,REF,ALT,INFO/AF,INFO/AF_nfe \
  /lustre/home/enza/franco/deliver/samples/napoli/embryos/embryos.chr21.gatk.vcf.gz \
  -Oz -o /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.vcf.gz

bcftools index /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.vcf.gz

(quick check)
bcftools view -h /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.vcf.gz | grep -E 'gnomAD_AF|gnomAD_AF_nfe'

STEP 3 — Re-run GREEN-VARAN with gnomAD present  [INSIDE CONTAINER, one-liner]
Adds: prioritization may activate (>0) if config thresholds met
apptainer exec -B /lustre/home/enza/franco/green-varan/resources:/opt/green-varan/resources -B /lustre/home/enza/franco/deliver:/deliver /lustre/home/enza/franco/green-varan/greenvaran_1.3.3.sif greenvaran smallvars --invcf /deliver/noncod_prior/embryos.chr21.gnomad.vcf.gz --outvcf /deliver/noncod_prior/embryos.chr21.greenvaran.vcf.gz --db /opt/green-varan/resources/GRCh38/GRCh38_GREEN-DB.bed.gz --dbschema /opt/green-varan/config/greendb_schema_v2.5.json --config /opt/green-varan/config/prioritize_smallvars.json

STEP 4 — Add ReMM score  [HOST, one-liner]
Adds: INFO/ReMM
# Añade la definición del tag y anota desde el TSV (3ra col = ReMM)
bcftools annotate \
  -a /lustre/home/enza/franco/green-varan/resources/GRCh38/GRCh38_ReMM.tsv.gz \
  -h <(echo '##INFO=<ID=ReMM,Number=1,Type=Float,Description="ReMM score (GRCh38)">') \
  -c CHROM,POS,ReMM \
  /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.vcf.gz \
  -Oz -o /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.remm.vcf.gz && \
bcftools index -f /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.remm.vcf.gz


STEP 5 — Re-run GREEN-VARAN with ReMM  [INSIDE CONTAINER, one-liner]
Adds: levels updated using AF + GREEN-DB + ReMM
apptainer exec -B /lustre/home/enza/franco/green-varan/resources:/opt/green-varan/resources -B /lustre/home/enza/franco/deliver:/deliver /lustre/home/enza/franco/green-varan/greenvaran_1.3.3.sif greenvaran smallvars --invcf /deliver/noncod_prior/embryos.chr21.gnomad.remm.vcf.gz --outvcf /deliver/noncod_prior/embryos.chr21.greenvaran.vcf.gz --db /opt/green-varan/resources/GRCh38/GRCh38_GREEN-DB.bed.gz --dbschema /opt/green-varan/config/greendb_schema_v2.5.json --config /opt/green-varan/config/prioritize_smallvars.json

STEP 6 — Add ncER percentile (BED)  [HOST, TWO-LINER]
Adds: INFO/ncER (noncoding intolerance)
echo '##INFO=<ID=ncER,Number=1,Type=Float,Description="ncER percentile (noncoding intolerance)">' > /lustre/home/enza/franco/deliver/noncod_prior/hdr.ncER.txt
bcftools annotate -h /lustre/home/enza/franco/deliver/noncod_prior/hdr.ncER.txt -a /lustre/home/enza/franco/green-varan/resources/GRCh38/GRCh38_ncER_perc.bed.gz -c CHROM,FROM,TO,INFO/ncER /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.remm.vcf.gz -Oz -o /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.remm.ncer.vcf.gz && bcftools index /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.remm.ncer.vcf.gz

STEP 7 — Re-run GREEN-VARAN with ncER  [INSIDE CONTAINER, one-liner]
Adds: levels updated using AF + ReMM + ncER
apptainer exec -B /lustre/home/enza/franco/green-varan/resources:/opt/green-varan/resources -B /lustre/home/enza/franco/deliver:/deliver /lustre/home/enza/franco/green-varan/greenvaran_1.3.3.sif greenvaran smallvars --invcf /deliver/noncod_prior/embryos.chr21.gnomad.remm.ncer.vcf.gz --outvcf /deliver/noncod_prior/embryos.chr21.greenvaran.vcf.gz --db /opt/green-varan/resources/GRCh38/GRCh38_GREEN-DB.bed.gz --dbschema /opt/green-varan/config/greendb_schema_v2.5.json --config /opt/green-varan/config/prioritize_smallvars.json

STEP 8 — Add FATHMM-MKL_NC (TSV)  [HOST, one-liner]
Adds: INFO/FATHMM_MKLNC
bcftools annotate \
  -a /lustre/home/enza/franco/green-varan/resources/GRCh38/GRCh38_FATHMM-MKL_NC.tsv.gz \
  -h <(echo '##INFO=<ID=FATHMM_NC,Number=1,Type=Float,Description="FATHMM-MKL non-coding score (GRCh38)">') \
  -c CHROM,POS,REF,ALT,FATHMM_NC \
  /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.gnomad.remm.ncer.vcf.gz \
  -Oz -o /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.fullannot.vcf.gz && \
bcftools index -f /lustre/home/enza/franco/deliver/noncod_prior/embryos.chr21.fullannot.vcf.gz


STEP 9 — Final GREEN-VARAN run (AF + ReMM + ncER + FATHMM)  [INSIDE CONTAINER, one-liner]
Adds: final prioritization (Levels L1–L4) using GREEN-DB + rarity + classic scores
apptainer exec -B /lustre/home/enza/franco/green-varan/resources:/opt/green-varan/resources -B /lustre/home/enza/franco/deliver:/deliver /lustre/home/enza/franco/green-varan/greenvaran_1.3.3.sif greenvaran smallvars --invcf /deliver/noncod_prior/embryos.chr21.fullannot.vcf.gz --outvcf /deliver/noncod_prior/embryos.chr21.greenvaran.final.vcf.gz --db /opt/green-varan/resources/GRCh38/GRCh38_GREEN-DB.bed.gz --dbschema /opt/green-varan/config/greendb_schema_v2.5.json --config /opt/green-varan/config/prioritize_smallvars.json

EXPECTED OUTPUTS
- embryos.chr21.gnomad.vcf.gz
- embryos.chr21.gnomad.remm.vcf.gz
- embryos.chr21.gnomad.remm.ncer.vcf.gz
- embryos.chr21.fullannot.vcf.gz
- embryos.chr21.greenvaran.final.vcf.gz (+ .tbi) and accompanying TSV/logs

QUICK CHECKS
- After each annotate: bcftools view -h <vcf.gz> | grep -E 'gnomAD_AF|ReMM|ncER|FATHMM_MKLNC'
- If see contig mismatches (chr21 vs 21), rename contigs on the sample VCF first, then repeat STEP 2.
 
