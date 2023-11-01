# squish

A program to re-arrange FASTQ file (.fastq.gz) reads in order to compress the file to a smaller size.

# Usage

Just point `squish` at your input fastq.gz file, and supply it with the filename of the desired output file.

```
$ ./squish data/test1.fastq.gz output.gz
2023/11/01 19:18:05 Input file data/test1.fastq.gz of size 4481 Bytes
2023/11/01 19:18:05 3276 reads loaded
2023/11/01 19:18:05 3276 reads sorted
2023/11/01 19:18:05 Output file created output.gz of size 4045 Bytes
2023/11/01 19:18:05 Size reduced by 436 Bytes (0.0973)
```

The current implementation of `squish` loads all fastq reads into memory, so you will need available system memory equal to the uncompressed size, plus extra for overhead. Hopefully future versions of `squish` can reduce memory requirements.

Tests with various .fastq.gz files saw compression improvements ranging from 3-35%.

# Installation

You can download a pre-built binary from the Releases page [here](https://github.com/stevekm/squish/releases).

You can build it from source with Go version 1.20+

```
go build -o ./squish ./main.go
```

(Makefile recipe for building is also included)



