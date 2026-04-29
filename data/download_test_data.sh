#!/bin/bash
set -eux

# download more fastq from here

# https://github.com/nf-core/test-datasets/tree/modules/data/genomics/homo_sapiens/illumina/fastq
# wget https://github.com/nf-core/test-datasets/blob/modules/data/genomics/homo_sapiens/illumina/fastq/test_1.fastq.gz

# wget https://github.com/nf-core/test-datasets/blob/modules/data/genomics/homo_sapiens/illumina/fastq/test2_2.fastq.gz

# wget https://github.com/nf-core/test-datasets/blob/modules/data/genomics/homo_sapiens/illumina/fastq/test_rnaseq_1.fastq.gz

# wget https://github.com/nf-core/test-datasets/blob/modules/data/genomics/homo_sapiens/illumina/fastq/test_rnaseq_2.fastq.gz





# SMALLER SIZE DATA SETS
# https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/samplesheet/v3.10/samplesheet_test.csv
# RAP1_IAA_30M_REP1,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_2.fastq.gz,reverse

for i in https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_1.fastq.gz https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_2.fastq.gz ; do wget $i ; done


# FULL SIZE DATA SETS
# https://github.com/nf-core/rnaseq/blob/master/conf/test_full.config
# https://raw.githubusercontent.com/nf-core/test-datasets/626c8fab639062eade4b10747e919341cbf9b41a/samplesheet/v3.10/samplesheet_full.csv
# s3://ngi-igenomes/test-data/rnaseq/SRX1603629_T1_1.fastq.gz,s3://ngi-igenomes/test-data/rnaseq/SRX1603629_T1_2.fastq.gz

for i in s3://ngi-igenomes/test-data/rnaseq/SRX1603629_T1_1.fastq.gz s3://ngi-igenomes/test-data/rnaseq/SRX1603629_T1_2.fastq.gz ; do aws s3 cp --no-sign-request $i ./ ; done


