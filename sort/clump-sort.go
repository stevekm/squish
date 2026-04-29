package sort

import (
	"bytes"
	go_sort "sort"
	fastq "squish/fastq"
)

type clumpRead struct {
	read fastq.FastqRead
	key  []byte
}

func SortReadsClump(reads *[]fastq.FastqRead) {
	SortReadsClumpK(reads, DefaultClumpKmerLen)
}

func SortReadsClumpK(reads *[]fastq.FastqRead, k int) {
	clumpReads := make([]clumpRead, len(*reads))
	for i, read := range *reads {
		clumpReads[i] = clumpRead{
			read: read,
			key:  ClumpKeyK(read, k),
		}
	}

	go_sort.Slice(clumpReads, func(i, j int) bool {
		return ClumpReadLess(clumpReads[i], clumpReads[j])
	})

	for i, clumpRead := range clumpReads {
		(*reads)[i] = clumpRead.read
	}
}

const DefaultClumpKmerLen = 16

func ClumpKey(read fastq.FastqRead) []byte {
	return ClumpKeyK(read, DefaultClumpKmerLen)
}

func ClumpKeyK(read fastq.FastqRead, k int) []byte {
	return clumpMinimizer(read.Sequence(), k)
}

func clumpMinimizer(sequence []byte, k int) []byte {
	sequence = trimLineEnding(sequence)
	if k < 1 {
		k = 1
	}
	if len(sequence) <= k {
		return canonicalKmer(sequence)
	}

	best := canonicalKmer(sequence[:k])
	for i := 1; i+k <= len(sequence); i++ {
		candidate := canonicalKmer(sequence[i : i+k])
		if bytes.Compare(candidate, best) < 0 {
			best = candidate
		}
	}
	return best
}

func canonicalKmer(kmer []byte) []byte {
	forward := append([]byte(nil), kmer...)
	reverseComplement := reverseComplement(kmer)
	if bytes.Compare(reverseComplement, forward) < 0 {
		return reverseComplement
	}
	return forward
}

func reverseComplement(sequence []byte) []byte {
	reverseComplement := make([]byte, len(sequence))
	for i, base := range sequence {
		reverseComplement[len(sequence)-1-i] = complementBase(base)
	}
	return reverseComplement
}

func complementBase(base byte) byte {
	switch base {
	case 'A', 'a':
		return 'T'
	case 'C', 'c':
		return 'G'
	case 'G', 'g':
		return 'C'
	case 'T', 't':
		return 'A'
	default:
		return 'N'
	}
}

func trimLineEnding(line []byte) []byte {
	line = bytes.TrimSuffix(line, []byte{'\n'})
	line = bytes.TrimSuffix(line, []byte{'\r'})
	return line
}

func ClumpCompare(a, b fastq.FastqRead) bool {
	return ClumpCompareK(a, b, DefaultClumpKmerLen)
}

func ClumpCompareK(a, b fastq.FastqRead, k int) bool {
	return ClumpReadLess(
		clumpRead{read: a, key: ClumpKeyK(a, k)},
		clumpRead{read: b, key: ClumpKeyK(b, k)},
	)
}

func ClumpReadLess(a, b clumpRead) bool {
	// Group reads by a canonical minimizer k-mer. This is closer to the
	// Clumpify-style goal of putting reads with shared sequence content near
	// each other, even when the shared sequence is not at the start of the read.
	if c := bytes.Compare(a.key, b.key); c != 0 {
		return c < 0
	}
	// Within a clump, keep records deterministic and compression-friendly by
	// sorting on the full sequence and quality bytes.
	if c := bytes.Compare(a.read.Sequence(), b.read.Sequence()); c != 0 {
		return c < 0
	}
	if c := bytes.Compare(a.read.QualityScores(), b.read.QualityScores()); c != 0 {
		return c < 0
	}
	if c := bytes.Compare(a.read.Id(), b.read.Id()); c != 0 {
		return c < 0
	}
	return a.read.I < b.read.I
}
