# Generating Whole-Genome LDZip Matrix from UK Biobank LD Matrices

This tutorial walks through generating whole-genome linkage disequilibrium (LD) matrices from pre-computed UK Biobank LD matrices (UKBB-LD) using the Nextflow pipeline provided in this repository.

**See also:** [1000 Genomes tutorial](g1k-tutorial.md) for generating LD matrices from raw VCF files.

---

## Setup: Clone and Build

Before running the pipeline, clone the repository and build the required binaries:

```bash
# Clone LDZip repository
git clone git@github.com:23andMe/LDZip.git
cd LDZip

# Build C++ binary
cd cpp
make
cd ..

# Install R package (installs dependencies, builds, and installs)
cd R
make install
cd ..

# Return to working directory
cd ..
```

After setup, the `ldzip` binary will be at `LDZip/cpp/bin/ldzip` and the R package `LDZipMatrix` will be installed.

---

## Quick Start: 10-Minute Test Run

Before downloading all 2,763 files, test with the first 4 chromosome 1 files:

```bash
mkdir -p data/ukbb

# Download first 4 files for quick test (chr1 chunks: 0-12Mb)
S3_BUCKET="s3://broad-alkesgroup-ukbb-ld/UKBB_LD"
chunks=(chr1_1_3000001 chr1_3000001_6000001 chr1_6000001_9000001 chr1_9000001_12000001)

for chunk in "${chunks[@]}"; do
  aws s3 cp --no-sign-request "${S3_BUCKET}/${chunk}.npz" data/ukbb/ &
  aws s3 cp --no-sign-request "${S3_BUCKET}/${chunk}.gz" data/ukbb/ &
done

wait

# Create test config (ukbb_test.yaml)
cat > ukbb_test.yaml << 'EOF'
npz_template: '${launchDir}/data/ukbb/chr{CHR}_{CHUNK}.npz'
outdir: 'output_test'
prefix: 'european_ukbb_test'
chroms: '1'
npz_ld_type: 'UNPHASED_R'
concat_pairwise: true
min: 0.1
bits: 8
EOF

# Run test pipeline
nextflow run LDZip/pipelines/wholeGenomeLD/main.nf -params-file ukbb_test.yaml -resume
```

This completes in ~10 minutes and produces a compressed LD matrix for the first 12Mb of chr1 in `output_test/whole_genome/`.

---

## Full Pipeline Overview

UKBB-LD provides summary linkage disequilibrium (LD) matrices computed from UK Biobank based on N=337K British-ancestry individuals. The LD matrices were computed by the [Alkes Price group at Harvard](https://labs.icahn.mssm.edu/minervalab/resources/data-ark/ukbb_ld/) and are publicly available in AWS S3. The LD information is stored as 2,763 3Mb-long regions spanning the entire genome in NPZ format.

The workflow consists of two main steps:

1. Download pre-computed NPZ LD matrices from AWS S3
2. Run a Nextflow pipeline to compress and concatenate into a whole genome LDZip matrix

---

## Prerequisites

- [`nextflow`](https://www.nextflow.io/docs/latest/install.html)
- [`ldzip`](../../README.md#installation) (C++ binary from this repository)
- AWS CLI (for downloading data)
- Sufficient compute (multi-core recommended)
- Adequate storage for final/intermediate files (~15GB for final output)

---

## Step 1: Download UKBB-LD NPZ Files

```bash
mkdir -p data/ukbb

# Download all NPZ files from AWS S3 (no AWS credentials required)
aws s3 sync --no-sign-request \
  s3://broad-alkesgroup-ukbb-ld/UKBB_LD/ \
  data/ukbb/
```

This will download 2,763 NPZ files (one per 3Mb region) with accompanying `.gz` variant metadata files.

---

## Step 2: Run Nextflow LD Pipeline

- Clone the LDZip repository:

  ```bash
  git clone git@github.com:23andMe/LDZip.git
  ```

- Create YAML file of input parameters

  **File: ukbb.yaml**
  ```yaml
  npz_template: '${launchDir}/data/ukbb/chr{CHR}_{CHUNK}.npz'
  outdir: 'output'
  prefix: 'european_ukbb'
  chroms: '20,21,22'
  npz_ld_type: 'UNPHASED_R'
  concat_pairwise: true
  min: 0.1
  bits: 8
  ```

  **Note:** The `npz_template` should point to where you downloaded the files in Step 1. Adjust the path as needed.

- If `ldzip` is NOT available in your `$PATH`, export its path before running the pipeline:

  ```bash
  export LDZIP=$(pwd)/LDZip/cpp/bin/ldzip
  ```

- Run the pipeline:

  ```bash
  nextflow run LDZip/pipelines/wholeGenomeLD/main.nf -params-file ukbb.yaml -resume
  ```

  **Tip:** For detailed parameter descriptions and validation rules, run:
  ```bash
  nextflow run LDZip/pipelines/wholeGenomeLD/main.nf --help
  ```

---

- **Note**

  The above configuration runs locally for chromosomes 20-22 only (`chroms`) so the pipeline completes quickly (typically ~10 minutes). This serves as a sanity check to ensure the inputs are correct and the workflow completes end-to-end with reasonable outputs.

  To run LDZip on all chromosomes, update to:

  ```yaml
  chroms: '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'
  ```

  **Important:** Set `concat_pairwise: true` for this dataset because the 3Mb regions have overlapping boundaries that cannot be handled by direct concatenation.

- **Running on HPC**

  For whole-genome processing, you may want to use an HPC cluster. For example, if using SLURM:

  **File: slurm.config**
  ```groovy
  params.partition        = "example_partition"

  process.executor        = "slurm"
  process.cpus            = 1
  process.errorStrategy   = 'retry'
  process.maxRetries      = 5
  process.queue           = params.partition

  executor.perCpuMemAllocation = true
  ```

  ```bash
  nextflow run LDZip/pipelines/wholeGenomeLD/main.nf -params-file ukbb.yaml -C slurm.config -resume
  ```

  For other HPC environments, refer to the Nextflow executor [guidelines](https://www.nextflow.io/docs/latest/executor.html).

---

## Expected Output

After successful completion, the final output will contain:

```
output/whole_genome/
├── european_ukbb.i.bin              # Row indices (compressed)
├── european_ukbb.i.bin.index        # Index for row indices
├── european_ukbb.meta.json          # Metadata
├── european_ukbb.p.bin              # Variant positions
├── european_ukbb.sqlite             # Variant index database
├── european_ukbb.vars.txt           # Variant information
├── european_ukbb.x.UNPHASED_R.bin   # LD values (compressed)
└── european_ukbb.x.UNPHASED_R.bin.index  # Index for LD values
```

**File sizes for whole genome (all chromosomes):**
```
5.6G    european_ukbb.i.bin
912K    european_ukbb.i.bin.index
4.0K    european_ukbb.meta.json
148M    european_ukbb.p.bin
1.2G    european_ukbb.sqlite
510M    european_ukbb.vars.txt
5.6G    european_ukbb.x.UNPHASED_R.bin
908K    european_ukbb.x.UNPHASED_R.bin.index
---
13G     total
```

---

## Using the LDZip Matrix in R

Once the pipeline completes, you can query the compressed LD matrix using the R package:

```r
library(LDZipMatrix)

# Load the LD matrix
ld <- LDZipMatrix("output/whole_genome/european_ukbb")

# Query LD between two variants
fetchLD(ld, "rs133036", "rs6001980")
# [1] 0.3700787

# Benchmark query time
system.time(fetchLD(ld, "rs133036", "rs6001980"))
#    user  system elapsed
#   0.003   0.000   0.010
```

---

## Parameters

| Parameter        | Description                                                                 | Example / Notes |
|------------------|-----------------------------------------------------------------------------|-----------------|
| `npz_template`   | Template path to per-chromosome-chunk NPZ files. `{CHR}` and `{CHUNK}` are replaced at runtime | `${launchDir}/data/ukbb/chr{CHR}_{CHUNK}.npz` |
| `outdir`         | Directory where LD outputs will be written                                | `output` |
| `prefix`         | Prefix used for naming output files                                       | `european_ukbb` |
| `chroms`         | Comma-separated list of chromosomes to process                            | `1–22` |
| `npz_ld_type`    | LD statistic type to compress from NPZ files                              | `UNPHASED_R`, `PHASED_R`, `UNPHASED_R2`, `DPRIME` |
| `concat_pairwise`| Use pairwise concatenation for overlapping regions (required for UKBB-LD) | `true` |
| `min`            | Minimum absolute LD value to store (values below are set to zero)        | `0.1` |
| `bits`           | Number of bits for quantization (8, 16, 32, or 99 for raw float)         | `8` |
| `stage_chunk`    | *(Optional)* Stage intermediate chunk files to outdir (default: `false`)  | `true` |
| `stage_chr`      | *(Optional)* Stage per-chromosome files to outdir (default: `false`)      | `true` |
