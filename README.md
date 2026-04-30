# squish

`squish` reorders FASTQ records before gzip compression. The goal is to place
similar records near each other so the compressed `.fastq.gz` output is smaller.

The current implementation supports in-memory sorting and external file-backed
bucket sorting. External sorting is the CLI default because real FASTQ files are
often larger than available RAM.

## Features

- Sort FASTQ records by sequence, GC content, quality string, heap-based
  sequence order, or clump-style minimizer grouping.
- External bucket sorting for large inputs.
- Arena-backed FASTQ parsing to reduce per-read allocation overhead.
- Optional paired/companion FASTQ reordering: sort R1 once, then apply the same
  `order.txt` permutation to R2 or other companion FASTQs.
- Optional mate-name validation for paired FASTQs.
- JSON run reports with read counts, file sizes, compression ratios, paths,
  bucket details, and paired output details.
- CPU and memory profiles for each run.
- Nextflow workflow for running multiple samples and sorting methods.
- Dockerfile for building a linux/amd64 container suitable for batch systems.

## Installation

Build locally with Go:

```bash
make build
```

or:

```bash
go build -trimpath -o ./squish ./cmd/squish
```

Run tests:

```bash
make test
```

## Go Package

The reusable implementation lives in the root Go package:

```go
import "squish"
```

External Go code can call `squish.Run` directly instead of shelling out to the
CLI:

```go
result, err := squish.Run(ctx, squish.Config{
    SortMethod:        "clump",
    SortEngine:        "external",
    InputFilepath:     "data/sample_R1.fastq.gz",
    OutputFilenameArg: "sample_R1.clump.fastq.gz",
    OutputDir:         "output",
    BucketStrategy:    "auto",
    BucketCount:       512,
    ClumpKmerLen:      31,
    CheckPairs:        true,
})
```

The CLI in `cmd/squish` is intentionally thin: it parses flags, builds a
`squish.Config`, calls `squish.Run`, and handles process exit codes.

## CLI Usage

Basic usage:

```bash
./squish [options] input.fastq.gz output.fastq.gz
```

By default, output files are written under `-outdir output`, so this command:

```bash
./squish data/test1.fastq.gz sorted.test1.fastq.gz
```

writes:

```text
output/sorted.test1.fastq.gz
output/order.txt
output/report.json
output/manifest.txt
output/profile.alpha/cpu.prof
output/profile.alpha/mem.prof
```

External sorting is the default:

```bash
./squish \
  -engine external \
  -m clump \
  -bucket auto \
  -buckets 512 \
  -clumpK 31 \
  -outdir output \
  data/sample_R1.fastq.gz \
  sample_R1.clump.fastq.gz
```

Use memory mode for smaller files or debugging:

```bash
./squish -engine memory -m alpha data/sample.fastq.gz sample.alpha.fastq.gz
```

## Sorting Methods

Use `-m` to choose the method:

- `alpha`: bytewise sequence sort.
- `gc`: sort by GC content.
- `qual`: sort by quality string.
- `clump`: groups reads by a pivot k-mer for compression-oriented clustering.
  For each read, every k-mer window is canonicalized (min of forward and
  reverse complement) and hashed with FNV-1a. The window with the highest
  hash value becomes the sort key. Reads sharing that k-mer land in the same
  clump. Using the max-hash pivot (rather than the lex-minimum) avoids
  over-grouping on low-complexity k-mers such as poly-A runs.

Tune the k-mer length with `-clumpK`. Longer k-mers are more specific and
give tighter clumps:

```bash
-clumpK 31
```

The default is `31`.

## External Buckets

Use `-bucket` to choose the external bucket strategy:

- `auto`: choose the default bucket strategy for the selected sorter.
- `sequence-prefix`: ordered buckets by sequence prefix.
- `quality-prefix`: ordered buckets by quality prefix.
- `gc-range`: ordered buckets by GC range.
- `hash`: fixed-count hash buckets.
- `clump-minimizer`: hash buckets based on the clump minimizer key.

`-buckets` controls the bucket count for configurable strategies:

```bash
./squish -engine external -bucket gc-range -buckets 1024 ...
```

Temporary bucket files are created under the output directory, in `tmp/<method>`
by default.

## Paired FASTQ Inputs

For paired-end data, R1 should define the sort order. `squish` writes `order.txt`
for the sorted R1 output, then companion FASTQs can be reordered using that same
permutation.

Example:

```bash
./squish \
  -engine external \
  -m clump \
  -paired data/sample_R2.fastq.gz \
  -pairedOut sample_R2.clump.fastq.gz \
  data/sample_R1.fastq.gz \
  sample_R1.clump.fastq.gz
```

Multiple companion files can be passed with commas or semicolons:

```bash
-paired data/sample_R2.fastq.gz,data/sample_I1.fastq.gz \
-pairedOut sample_R2.sorted.fastq.gz,sample_I1.sorted.fastq.gz
```

Pair checking is enabled by default:

```bash
-checkPairs=true
```

The check normalizes read IDs by removing `@`, whitespace suffixes, and trailing
`/1` or `/2`. Disable it only when read names are known to differ but record
order is guaranteed to match:

```bash
-checkPairs=false
```

Companion FASTQ reordering uses a temporary record file plus offsets, so it does
not keep the full companion FASTQ in memory.

## Report Output

Every run writes a JSON report, default `report.json`, under the output
directory. It includes:

- version, start/end time, and duration
- sort method, engine, bucket strategy, and clump k-mer length
- input and output file paths and sizes
- read counts and uncompressed bytes processed
- output compression ratio and size reduction ratio
- profile paths
- manifest path
- paired/companion FASTQ output details when `-paired` is used

`manifest.txt` is a simple newline-delimited list of absolute paths to every
output FASTQ produced by the run. It includes the primary sorted FASTQ and any
paired/companion FASTQs.

Example report fields:

```json
{
  "sort_method": "clump",
  "sort_engine": "external",
  "clump_kmer_length": 31,
  "reads": 1000000,
  "output": {
    "path": "/absolute/path/output/sample_R1.clump.fastq.gz",
    "argument": "sample_R1.clump.fastq.gz"
  },
  "manifest_file": {
    "path": "output/manifest.txt"
  },
  "paired_outputs": [
    {
      "input": { "path": "data/sample_R2.fastq.gz" },
      "output": { "argument": "sample_R2.clump.fastq.gz" }
    }
  ]
}
```

## Nextflow Workflow

The workflow in `main.nf` runs every sample in a samplesheet through every method
listed in `params.methods`.

Default settings live in `nextflow.config`. The default samplesheet is currently:

```groovy
params.samplesheet = 'samples.small.csv'
```

Run:

```bash
nextflow run main.nf
```

or override inputs:

```bash
nextflow run main.nf \
  --samplesheet samples.csv \
  --engine external \
  --methods alpha,gc,qual,clump \
  --clump_k 31
```

The workflow expects the `squish` binary to be available on `PATH` inside the
task environment. When using Docker, the container should include `squish`.

### Samplesheet Format

Samplesheets use this header:

```csv
sample_id,fastqin,fastqout,paired_fastqs,paired_outputs
```

Columns:

- `sample_id`: output grouping name.
- `fastqin`: primary FASTQ, usually R1.
- `fastqout`: output filename base for the primary sorted FASTQ.
- `paired_fastqs`: optional semicolon-separated companion FASTQs.
- `paired_outputs`: optional semicolon-separated output filename bases for
  companion FASTQs.

Example:

```csv
sample_id,fastqin,fastqout,paired_fastqs,paired_outputs
SRR6357076,data/SRR6357076_1.fastq.gz,squish.SRR6357076_1,data/SRR6357076_2.fastq.gz,squish.SRR6357076_2
```

For each method, outputs are published under:

```text
output-nextflow/<sample_id>/<method>/
```

The workflow also writes:

- `output-nextflow/report.csv`
- `output-nextflow/nextflow-report.html`
- `output-nextflow/nextflow-trace.txt`
- `output-nextflow/nextflow-timeline.html`

Each sample/method output directory also contains `manifest.<method>.txt`.

## Docker

Build a Docker image:

```bash
make docker-build
```

The Dockerfile is pinned for linux/amd64 builds, which is useful when building
on Apple Silicon but running in x86 Linux environments. Nextflow itself is not
installed in the container; Nextflow runs the container.

## Test Data

The `data/download_test_data.sh` script downloads the small and full-size test
FASTQs used during development and recreates derived subset files.

```bash
bash data/download_test_data.sh
```

The script skips files that already exist.

## Development

Common commands:

```bash
make test
make build
make test-run-all
make test-run-all-external
make nextflow-test-run-all-external
```

`make test-run-all` runs all sort methods on the configured `FASTQIN`.

# AI Usage

Releases after v0.2 (c385a1) were developed with Codex 5.5.

# References

Inspired by the concepts behind Clumpify by Brian Bushnell:

- https://www.biostars.org/p/225338/
- https://github.com/BioInfoTools/BBMap/blob/master/current/clump/Clumpify.java
- https://github.com/BioInfoTools/BBMap/blob/master/sh/clumpify.sh

Inspired by RustQC (https://seqera.io/blog/rustqc/) and https://rewrites.bio/ by Phil Ewels.
