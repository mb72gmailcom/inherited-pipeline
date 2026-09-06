# Inherited

Analyze rare variants in family trios from a VCF, precomputed gnomAD allele-frequency JSON, and a family relations file.

## Run

**Without installing** (clone and run from the project directory):

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22
```

**Or install** so the `inherited` command is on your PATH:

```bash
pip install -e ".[dev]"
inherited analyze --vcf chr22.vcf.gz --af-json gnomad_chr22.json --family-file families.tsv -o results/chr22
```

**Split pVCFs** (several `{start}_{end}` files per chromosome):

```bash
python run.py analyze \
  --vcf-dir /data/spark/vcf/chr2 \
  --vcf-pattern SPARK.WGS.2026_08.gatk \
  --af-json gnomad_chr2.json \
  --family-file families.tsv \
  -o results/chr2
```

`--vcf-dir` looks for `{pattern}.{chr}_{start}_{end}.vcf.gz` (or `.vcf`) in a directory that holds one chromosome. `--vcf-pattern` may be the callset stem (`SPARK.WGS.2026_08.gatk`) or that stem plus the contig (`...gatk.chr21`). Output files are labeled by each shard's coordinates (`inherited_95000001_97500000.tsv`), not by `--segment-size`. `--vcf` and `--vcf-dir` are mutually exclusive.

To analyze one already-split shard with `--vcf` and write labeled files into a shared output directory:

```bash
python run.py analyze \
  --vcf SPARK.WGS.2026_08.gatk.chr22_95000001_97500000.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22 \
  --file-suffix 95000001_97500000
```

`--file-suffix` requires `--vcf`, turns off output segmentation, and appends `_{suffix}` to every output file (`inherited_95000001_97500000.tsv`, `stats_95000001_97500000.json`, `params_95000001_97500000.json`, and so on). It cannot be used with `--vcf-dir` or `--resume`.

## Family file format

Tab-separated file with a header row (required via `--family-file`):

```text
spid    sfid    father    mother    sex
child1  fam1    fa1       ma1       Female
boy1    fam1    fa1       ma1       1
fa1     fam1    0         0         Male
ma1     fam1    0         0         2
```

Required fields (either name is accepted): `spid` / `ind_id`, `sfid` / `family_id`, `father` / `father_id`, `mother` / `mother_id`, `sex`. Extra columns are ignored.

If a `sample_id` column is present, child / father / mother IDs used for VCF matching are remapped `ind_id → sample_id` from rows that have both values. Hits in the TSV and `*_per_person.json` files are labeled with `ind_id` (via the inverse map). A complete trio is dropped when any of the three has no mapping (for example a parent who was not sequenced). Files without `sample_id` keep IDs as written.

If a file uses other header names, pass `--family-map` with a JSON object from the internal name to that header:

```json
{"spid": "IID", "sfid": "FID", "father": "PAT", "mother": "MAT", "sex": "sex"}
```

Unmapped fields still use the aliases above. If both names of an alias pair are present (e.g. `spid` and `ind_id`), `--family-map` is required for that field.

Complete trios (`father` and `mother` both present — not empty, `0`, or `False`) with a recognized sex are used for analysis. Sex values: males `1` / `Male` / `male`; females `2` / `Female` / `female`. Children with missing or unrecognized sex are skipped on all chromosomes.

## gnomAD common-af.json

```json
{
  "var_key": 0.1,
  "chr22_12345_A_G" : 0.1,...
}
```

Variants with AF above the threshold are skipped. Missing keys default to `0` and are kept.

## Usage

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22
```

Disable multiallelic per-ALT processing:

```bash
python run.py analyze --vcf chr22.vcf.gz --af-json gnomad.json --family-file families.tsv -o results/chr22 --no-multiallelic
```

Enable debug memory logging every 50,000 variants (custom interval with `--memory-block`):

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad.json \
  --family-file families.tsv \
  -o results/chr22 \
  --debug \
  --memory-block 50000
```

Skip variants inside repetitive regions (one file per chromosome):

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22 \
  --remove-repeats chr22_repeats.bed
```

Repeat files are 0-based BED-style ``chrom start end`` rows with half-open intervals ``[start, end)``. VCF positions are 1-based, so a row ``chr22 2500 3500`` skips POS 2501 through 3500.

## Quality filters

Diploid QC (autosomes, female chrX, male PAR): `DP≥10`, `GQ≥20`. Allele balance applies only when the genotype carries the queried alt: heterozygotes (`ac==1`) require `0.25 ≤ AB ≤ 0.75`; homozygous-alt (`ac≥2`) require `AB≥0.9`. Homozygous-reference calls skip AB unless `--ab-hom00-threshold` is set; then diploid `0/0` requires `ab00 ≥ value` where `ab00 = AD[REF] / sum(AD)` (typical `0.9`).

Haploid QC (male chrX nonPAR, chrY): `DP≥5`, `GQ≥20`, and `AB≥0.85` when the genotype carries the alt.

Override any of these from the command line (omitted flags use the defaults above):

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22 \
  --gq-threshold 20 \
  --dp-threshold 10 \
  --dp-haploid-threshold 5 \
  --ab-threshold 0.25 \
  --ab-hom-threshold 0.9 \
  --ab-haploid-threshold 0.85
```

`--ab-threshold` is the diploid heterozygous half-band: het AB must fall in `[value, 1-value]`. `--ab-hom00-threshold` is an optional diploid `0/0` REF-fraction floor (`ab00 ≥ value`); omit it to keep the historical no-AB-on-`0/0` behavior. Applied values are written to `params.json` under `quality_filters`.

## Chromosome X

VCFs are expected to be split per chromosome. Contigs `X` and `chrX` use sex-aware logic:

- Female children: same diploid QC and inheritance rules as autosomes; one output bucket
- Male PAR1 / PAR2 (hg38 closed intervals `[10001, 2781479]` and `[155701383, 156030895]`): same diploid trio logic as autosomes
- Male nonPAR: mother–son pairs only (father ignored on chrX); child uses haploid QC (`DP≥5`, `AB≥0.85`, `GQ≥20`); mother uses diploid QC; inherited if `ac>0` and `mac>0`; denovo if `ac>0` and `mac==0` with mother genotype homozygous reference

Diploid call classes (autosomes, female chrX, male PAR):

- Skip the trio if either parent fails QC (`mac < 0` or `fac < 0`)
- denovo: `ac>0` and `mac==0` and `fac==0`, both parents `0/0`, and the child genotype uses only ref and this alt (e.g. `0/1` or `1/1`). Child `1/2` from `0/0 × 0/0` is `mendelian_bad`, not denovo. Parents that lack this alt but carry another allele (e.g. `0/2` when testing allele 3) are skipped
- inherited: at least one parent carries this alt and the full trio genotypes are Mendelian-compatible; stored parental genotypes are the real calls (no synthetic `0/0`)
- mendelian_bad: the child genotype is not Mendelian given both parents (including `0/2 × 0/2 → 1/2`, and `0/0 × 0/0 → 1/2`)

## Chromosome Y

Contigs `Y` and `chrY` are treated as nonPAR only. Analysis loops over male children and father–son pairs:

- Child and father both use haploid QC (`DP≥5`, `AB≥0.85`, `GQ≥20`)
- Inherited if `ac>0` and `fac>0`; denovo if `ac>0` and `fac==0` with father genotype homozygous reference
- Father haploid QC failure skips the pair
- Output uses only the `males_nonPar` bucket (no empty PAR/female files)
- Full format: `child_id=father_gt|child_gt|child_gq`

## Output

Autosomal output directory contains:

- `inherited_XXXXX.tsv` / `mendelian_bad_XXXXX.tsv` / `denovo_XXXXX.tsv` — segmented result files (when `--vcf` and `--segment-size > 0`)
- `inherited_{start}_{end}.tsv` / `mendelian_bad_{start}_{end}.tsv` / `denovo_{start}_{end}.tsv` — one file per input shard (when `--vcf-dir`)
- `inherited.tsv` / `mendelian_bad.tsv` / `denovo.tsv` — single files when `--vcf` and `--segment-size 0`
- `inherited_{suffix}.tsv` / `mendelian_bad_{suffix}.tsv` / `denovo_{suffix}.tsv` — single files when `--vcf` and `--file-suffix` (every sidecar JSON is labeled the same way)
- `inherited_per_variant.json`, `inherited_per_person.json`, `denovo_per_variant.json`, `denovo_per_person.json`, `mendelian_bad_per_gt.json`, `stats.json`

chrX output uses sex/region buckets:

- `inherited_females_XXXXX.tsv`, `inherited_males_par1_XXXXX.tsv`, `inherited_males_nonPar_XXXXX.tsv`, `inherited_males_par2_XXXXX.tsv`
- matching `mendelian_bad_*` and `denovo_*` files
- matching per-bucket `inherited_per_person_*.json`, `denovo_per_person_*.json`, `mendelian_bad_per_gt_*.json`, `inherited_per_variant_*.json`, `denovo_per_variant_*.json`, `stats_*.json`

Shared run metadata:

- `checkpoint.json` — resume point (updated after each segment)
- `stats_cumulative.json` — running totals (updated after each segment)
- `cumulative_detail_deltas.jsonl` — append-only per-segment summary details used for resume
- `stats.json` — final overall summary counts
- `params.json` — parameters used for this run

Each chromosome output directory (e.g. `results/chr22/`) contains its own `params.json` alongside the result files.

Result TSV columns (default `--short-format`):

```text
#CHROM  POS  ID  REF  ALT  PATIENTS
22      3000 .   A   G    child1;child2
```

Use `--no-short-format` for full genotype output:

```text
#CHROM  POS  ID  REF  ALT  TRIO_CALLS
22      3000 .   A   G    child1=0/1|0/0|0/1|30
```

Male nonPAR chrX full format omits the father genotype: `child_id=mother_gt|child_gt|child_gq`.

Male nonPAR chrY full format omits the mother genotype: `child_id=father_gt|child_gt|child_gq`.

Use `--block-size` (default `10000`) for in-memory buffer flushes within a segment.

Use `--segment-size` (default `1000000`) with `--vcf` to split output into segment files. Set `--segment-size 0` to disable segmentation. `--segment-size` is ignored with `--vcf-dir`, where each input shard writes files labeled `{start}_{end}`. `--file-suffix` also disables segmentation.

Resume after a crash:

```bash
python run.py analyze ... -o results/chr2 --resume
```

Requires an existing incomplete `checkpoint.json`. With `--vcf`, also requires `--segment-size > 0`. With `--vcf-dir`, resume skips shards that ended at or before the checkpointed position.

## Tests

```bash
pytest
```
