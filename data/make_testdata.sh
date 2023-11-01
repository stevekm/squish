#!/bin/bash
# inflate a larger test data file
set -eu

for i in $(seq 1 200); do
    cat test1.fastq >> test2.fastq
done

gzip -9 -c test2.fastq > test2.fastq.gz