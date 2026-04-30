package sort

import (
	"bytes"
	go_sort "sort"
	fastq "squish/fastq"
)

type clumpRead struct {
	read fastq.FastqRead
	// key is precomputed once per read so sort comparisons do not rescan the
	// sequence on every Less call. This is important because sort comparisons
	// happen O(n log n) times.
	key []byte
}

func SortReadsClump(reads *[]fastq.FastqRead) {
	SortReadsClumpK(reads, DefaultClumpKmerLen)
}

func SortReadsClumpK(reads *[]fastq.FastqRead, k int) {
	// Build a sidecar slice containing each read plus its minimizer key. The
	// final loop copies the sorted read order back into the caller's slice.
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
	// Use the lexicographically smallest canonical k-mer as the clump key.
	// Canonicalization compares the forward k-mer to its reverse complement so
	// reads from opposite strands can still land near each other.
	if k < 1 {
		k = 1
	}
	if len(sequence) <= k {
		return canonicalKmer(sequence)
	}

	// Preallocate one reverse-complement buffer and one best-key buffer for the
	// entire read so the inner loop produces zero heap allocations.
	rcBuf := make([]byte, k)
	best := make([]byte, k)

	reverseComplementInto(sequence[:k], rcBuf)
	if bytes.Compare(rcBuf, sequence[:k]) < 0 {
		copy(best, rcBuf)
	} else {
		copy(best, sequence[:k])
	}

	for i := 1; i+k <= len(sequence); i++ {
		kmer := sequence[i : i+k]
		reverseComplementInto(kmer, rcBuf)
		var canonical []byte
		if bytes.Compare(rcBuf, kmer) < 0 {
			canonical = rcBuf
		} else {
			canonical = kmer
		}
		if bytes.Compare(canonical, best) < 0 {
			copy(best, canonical)
		}
	}
	return best
}

func canonicalKmer(kmer []byte) []byte {
	rc := reverseComplement(kmer)
	if bytes.Compare(rc, kmer) < 0 {
		return rc
	}
	return append([]byte(nil), kmer...)
}

func reverseComplement(sequence []byte) []byte {
	rc := make([]byte, len(sequence))
	reverseComplementInto(sequence, rc)
	return rc
}

func reverseComplementInto(src, dst []byte) {
	for i, base := range src {
		dst[len(src)-1-i] = complementBase(base)
	}
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
		// Ambiguous bases are deliberately retained as an ambiguous placeholder
		// instead of trying to force them into a two-bit A/C/G/T encoding.
		return 'N'
	}
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
