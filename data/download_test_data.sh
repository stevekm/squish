#!/bin/bash
set -euo pipefail

# Download and derive the FASTQ files used for local squish testing.
#
# This script is intentionally idempotent: existing files are left in place.
# Run it from anywhere; outputs are written next to this script in data/.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

download_url() {
    local url="$1"
    local output="$2"

    if [[ -s "${output}" ]]; then
        echo "SKIP existing ${output}"
        return
    fi

    echo "DOWNLOAD ${url} -> ${output}"
    if command -v curl >/dev/null 2>&1; then
        curl -L --fail -o "${output}" "${url}"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "${output}" "${url}"
    else
        echo "ERROR: curl or wget is required" >&2
        exit 1
    fi
}

download_s3() {
    local s3_path="$1"
    local output="$2"

    if [[ -s "${output}" ]]; then
        echo "SKIP existing ${output}"
        return
    fi
    if ! command -v aws >/dev/null 2>&1; then
        echo "ERROR: aws CLI is required to download ${s3_path}" >&2
        exit 1
    fi

    echo "DOWNLOAD ${s3_path} -> ${output}"
    aws s3 cp --no-sign-request "${s3_path}" "${output}"
}

make_head_fastq_gz() {
    local source_fastq="$1"
    local line_count="$2"
    local output_gz="$3"
    local output_fastq="${output_gz%.gz}"
    local tmp_gz="${output_gz}.tmp"

    if [[ -s "${output_gz}" ]]; then
        echo "SKIP existing ${output_gz}"
        return
    fi
    if [[ ! -s "${source_fastq}" ]]; then
        echo "ERROR: source FASTQ does not exist: ${source_fastq}" >&2
        exit 1
    fi

    echo "DERIVE head -n ${line_count} ${source_fastq} -> ${output_gz}"
    head -n "${line_count}" "${source_fastq}" > "${output_fastq}"
    gzip -c "${output_fastq}" > "${tmp_gz}"
    mv "${tmp_gz}" "${output_gz}"
    rm -f "${output_fastq}"
}

# Small nf-core module FASTQs currently present in data/.
download_url \
    "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_1.fastq.gz" \
    "test_1.fastq.gz"
download_url \
    "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test2_2.fastq.gz" \
    "test2_2.fastq.gz"
download_url \
    "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_rnaseq_1.fastq.gz" \
    "test_rnaseq_1.fastq.gz"
download_url \
    "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_rnaseq_2.fastq.gz" \
    "test_rnaseq_2.fastq.gz"

# Smaller paired RNA-seq test data from the nf-core/rnaseq test samplesheet.
# RAP1_IAA_30M_REP1:
# https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/samplesheet/v3.10/samplesheet_test.csv
download_url \
    "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_1.fastq.gz" \
    "SRR6357076_1.fastq.gz"
download_url \
    "https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_2.fastq.gz" \
    "SRR6357076_2.fastq.gz"


download_url \
"https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357073_1.fastq.gz" "SRR6357073_1.fastq.gz"

download_url \
"https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357074_1.fastq.gz" "SRR6357074_1.fastq.gz"

download_url \
"https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_1.fastq.gz" "SRR6357076_1.fastq.gz"

download_url \
"https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_2.fastq.gz" "SRR6357076_2.fastq.gz"

download_url \
"https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357075_1.fastq.gz" "SRR6357075_1.fastq.gz"

# Full-size paired RNA-seq test data from the nf-core/rnaseq full test config.
# https://github.com/nf-core/rnaseq/blob/master/conf/test_full.config
# https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/samplesheet/v3.10/samplesheet_full.csv
download_s3 \
    "s3://ngi-igenomes/test-data/rnaseq/SRX1603629_T1_1.fastq.gz" \
    "SRX1603629_T1_1.fastq.gz"
download_s3 \
    "s3://ngi-igenomes/test-data/rnaseq/SRX1603629_T1_2.fastq.gz" \
    "SRX1603629_T1_2.fastq.gz"

# Keep the decompressed read 1 FASTQ because the subset files below were made
# with head -n from the uncompressed file.
if [[ ! -s "SRX1603629_T1_1.fastq" ]]; then
    echo "DECOMPRESS SRX1603629_T1_1.fastq.gz -> SRX1603629_T1_1.fastq"
    gzip -dc "SRX1603629_T1_1.fastq.gz" > "SRX1603629_T1_1.fastq"
else
    echo "SKIP existing SRX1603629_T1_1.fastq"
fi

# Derived subsets currently present in data/. The number in the filename is
# the number of FASTQ text lines retained by head -n, not the number of reads.
make_head_fastq_gz "SRX1603629_T1_1.fastq" 400000 "SRX1603629_T1_1.400000.fastq.gz"
make_head_fastq_gz "SRX1603629_T1_1.fastq" 4000000 "SRX1603629_T1_1.4000000.fastq.gz"
make_head_fastq_gz "SRX1603629_T1_1.fastq" 40000000 "SRX1603629_T1_1.40000000.fastq.gz"

echo "Done."
