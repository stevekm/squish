package sort

import (
	"bytes"
	go_sort "sort"
	fastq "squish/fastq"
)

type clumpRead struct {
	read fastq.FastqRead
	// key is the canonical pivot k-mer, precomputed once so sort comparisons
	// do not rescan the sequence on every Less call (O(n log n) comparisons).
	key []byte
	// pivotPos is the position in the read where the pivot k-mer starts.
	// Within a clump (same key), sorting by ascending position places reads
	// whose shared k-mer begins at the same offset adjacent to each other,
	// giving LZ77 the longest possible run of identical bytes.
	pivotPos int
	// rcFlipped is true when the canonical pivot was the reverse complement of
	// the forward k-mer, meaning the read aligns to the minus strand. When
	// RComp is enabled, these reads are reverse-complemented in the output so
	// all reads in a clump are on the same strand.
	rcFlipped bool
}

// ClumpSortOptions controls the behaviour of SortReadsClumpOpts.
type ClumpSortOptions struct {
	K        int  // k-mer length; 0 falls back to DefaultClumpKmerLen
	MinCount int  // ignore pivot k-mers appearing fewer than MinCount times (0 = disabled)
	RComp    bool // reverse-complement reads whose pivot was on the minus strand
	RawPivot bool // pick the lex-max canonical k-mer instead of the max-hash k-mer
	Border   int  // number of bases excluded from each end of the read during pivot selection
}

// SortReadsClump sorts using default k and no extra options.
func SortReadsClump(reads *[]fastq.FastqRead) {
	SortReadsClumpOpts(reads, ClumpSortOptions{K: DefaultClumpKmerLen})
}

// SortReadsClumpK sorts with a custom k and no extra options.
func SortReadsClumpK(reads *[]fastq.FastqRead, k int) {
	SortReadsClumpOpts(reads, ClumpSortOptions{K: k})
}

// SortReadsClumpOpts is the full implementation that honours all options.
func SortReadsClumpOpts(reads *[]fastq.FastqRead, opts ClumpSortOptions) {
	k := opts.K
	if k < 1 {
		k = DefaultClumpKmerLen
	}

	// Build a k-mer frequency table when mincount filtering is requested.
	// Only k-mers appearing at least MinCount times are eligible as pivots,
	// which avoids grouping reads by error k-mers that occur only once.
	var eligible func([]byte) bool
	if opts.MinCount > 1 {
		counts := countKmers(*reads, k)
		minCount := opts.MinCount
		eligible = func(kmer []byte) bool {
			return counts[string(kmer)] >= minCount
		}
	}

	border := opts.Border
	if border < 0 {
		border = 0
	}

	clumpReads := make([]clumpRead, len(*reads))
	for i, read := range *reads {
		key, pos, rcFlipped := clumpMinimizerFull(read.Sequence(), k, opts.RawPivot, eligible, border)
		clumpReads[i] = clumpRead{
			read:      read,
			key:       key,
			pivotPos:  pos,
			rcFlipped: rcFlipped,
		}
	}

	go_sort.Slice(clumpReads, func(i, j int) bool {
		return ClumpReadLess(clumpReads[i], clumpReads[j])
	})

	for i, cr := range clumpReads {
		if opts.RComp && cr.rcFlipped {
			// Flip reads whose pivot was on the minus strand so all reads in a
			// clump are in the same orientation — consecutive sequence lines
			// become more byte-similar, increasing LZ77 back-reference density.
			cr.read.OverrideSeq = reverseComplement(cr.read.Sequence())
			cr.read.OverrideQual = reverseBytes(cr.read.QualityScores())
		}
		(*reads)[i] = cr.read
	}
}

const DefaultClumpKmerLen = 31
const DefaultClumpBorder = 1 // exclude outermost base on each end from pivot selection (matches Clumpify default)

func ClumpKey(read fastq.FastqRead) []byte {
	return ClumpKeyK(read, DefaultClumpKmerLen)
}

func ClumpKeyK(read fastq.FastqRead, k int) []byte {
	key, _, _ := clumpMinimizerFull(read.Sequence(), k, false, nil, 0)
	return key
}

// hashKmer returns a 64-bit FNV-1a hash of a k-mer byte slice. The hash is
// position-sensitive and well-distributed, which breaks the lex-minimum bias
// toward poly-A / low-complexity k-mers.
func hashKmer(kmer []byte) uint64 {
	const (
		fnvOffset = uint64(14695981039346656037)
		fnvPrime  = uint64(1099511628211)
	)
	h := fnvOffset
	for _, b := range kmer {
		h ^= uint64(b)
		h *= fnvPrime
	}
	return h
}

// clumpMinimizerFull returns the pivot k-mer, its position in the sequence,
// and whether the canonical form was the reverse complement (rcFlipped).
//
// When rawPivot is false (default), the pivot is chosen by max FNV-1a hash of
// the canonical k-mer. When rawPivot is true, the lex-maximum canonical k-mer
// is chosen instead — this clusters reads by nucleotide composition rather
// than by hash, which can slightly improve compression for some datasets.
//
// When eligible is non-nil, only k-mers for which eligible returns true are
// candidates. If no k-mer passes the filter, a nil key is returned (the read
// sorts into an unclustered group at the front).
func clumpMinimizerFull(sequence []byte, k int, rawPivot bool, eligible func([]byte) bool, border int) (key []byte, pos int, rcFlipped bool) {
	if k < 1 {
		k = 1
	}
	if border < 0 {
		border = 0
	}
	// If the border-trimmed region is too short, fall back to no border.
	if border > 0 && len(sequence)-2*border < k {
		border = 0
	}
	if len(sequence) <= k {
		canon := canonicalKmer(sequence)
		rc := reverseComplement(sequence)
		flipped := bytes.Compare(rc, sequence) < 0
		if eligible != nil && !eligible(canon) {
			return nil, 0, false
		}
		return canon, 0, flipped
	}

	rcBuf := make([]byte, k)
	best := make([]byte, k)
	bestPos := -1
	var bestHash uint64
	bestRC := false

	for i := border; i+k <= len(sequence)-border; i++ {
		kmer := sequence[i : i+k]
		reverseComplementInto(kmer, rcBuf)
		thisRC := bytes.Compare(rcBuf, kmer) < 0
		var canonical []byte
		if thisRC {
			canonical = rcBuf
		} else {
			canonical = kmer
		}
		if eligible != nil && !eligible(canonical) {
			continue
		}
		if bestPos == -1 {
			// First eligible k-mer initialises best.
			copy(best, canonical)
			bestPos = i
			bestRC = thisRC
			if !rawPivot {
				bestHash = hashKmer(best)
			}
			continue
		}
		var better bool
		if rawPivot {
			better = bytes.Compare(canonical, best) > 0
		} else {
			h := hashKmer(canonical)
			if better = h > bestHash; better {
				bestHash = h
			}
		}
		if better {
			copy(best, canonical)
			bestPos = i
			bestRC = thisRC
		}
	}

	if bestPos == -1 {
		// No eligible k-mer found — return nil key so the read sorts with
		// other unclustered reads rather than being given an arbitrary pivot.
		return nil, 0, false
	}
	return best, bestPos, bestRC
}

// countKmers builds a frequency table of all canonical k-mers across reads.
func countKmers(reads []fastq.FastqRead, k int) map[string]int {
	counts := make(map[string]int)
	rcBuf := make([]byte, k)
	for _, read := range reads {
		seq := read.Sequence()
		for i := 0; i+k <= len(seq); i++ {
			kmer := seq[i : i+k]
			reverseComplementInto(kmer, rcBuf)
			var canonical []byte
			if bytes.Compare(rcBuf, kmer) < 0 {
				canonical = rcBuf
			} else {
				canonical = kmer
			}
			counts[string(canonical)]++
		}
	}
	return counts
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

func reverseBytes(b []byte) []byte {
	out := make([]byte, len(b))
	for i, v := range b {
		out[len(b)-1-i] = v
	}
	return out
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
	keyA, posA, rcA := clumpMinimizerFull(a.Sequence(), k, false, nil, 0)
	keyB, posB, rcB := clumpMinimizerFull(b.Sequence(), k, false, nil, 0)
	return ClumpReadLess(
		clumpRead{read: a, key: keyA, pivotPos: posA, rcFlipped: rcA},
		clumpRead{read: b, key: keyB, pivotPos: posB, rcFlipped: rcB},
	)
}

func ClumpReadLess(a, b clumpRead) bool {
	// Primary: group reads by their pivot k-mer.
	if c := bytes.Compare(a.key, b.key); c != 0 {
		return c < 0
	}
	// Secondary: within a clump, sort by the position of the pivot k-mer.
	// Reads sharing the same k-mer at the same offset have identical bytes
	// at that position, giving LZ77 the longest possible back-references.
	if a.pivotPos != b.pivotPos {
		return a.pivotPos < b.pivotPos
	}
	// Tertiary: full sequence and quality bytes for deterministic output.
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
