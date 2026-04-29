package sort

import (
	"bytes"
	go_sort "sort"
	fastq "squish/fastq"
)

func SortReadsClump(reads *[]fastq.FastqRead) {
	go_sort.Slice(*reads, func(i, j int) bool {
		return ClumpCompare((*reads)[i], (*reads)[j])
	})
}

const clumpKmerLen = 16

func ClumpKey(read fastq.FastqRead) []byte {
	return clumpMinimizer(read.Sequence(), clumpKmerLen)
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
	// Group reads by a canonical minimizer k-mer. This is closer to the
	// Clumpify-style goal of putting reads with shared sequence content near
	// each other, even when the shared sequence is not at the start of the read.
	if c := bytes.Compare(ClumpKey(a), ClumpKey(b)); c != 0 {
		return c < 0
	}
	// Within a clump, keep records deterministic and compression-friendly by
	// sorting on the full sequence and quality bytes.
	if c := bytes.Compare(a.Sequence(), b.Sequence()); c != 0 {
		return c < 0
	}
	if c := bytes.Compare(a.QualityScores(), b.QualityScores()); c != 0 {
		return c < 0
	}
	if c := bytes.Compare(a.Id(), b.Id()); c != 0 {
		return c < 0
	}
	return a.I < b.I
}
