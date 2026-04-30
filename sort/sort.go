package sort

import (
	go_sort "sort"
	fastq "squish/fastq"
)

// using built-in Go sort methods
func SortReadsSequence(reads *[]fastq.FastqRead) {
	sorter := AlphaSort{}
	go_sort.Slice((*reads), func(i, j int) bool { return sorter.Less((*reads)[i], (*reads)[j]) })
}

func SortReadsGC(reads *[]fastq.FastqRead) {
	sorter := GCSort{}
	go_sort.Slice((*reads), func(i, j int) bool { return sorter.Less((*reads)[i], (*reads)[j]) })
}

func SortReadsQual(reads *[]fastq.FastqRead) {
	sorter := QualitySort{}
	go_sort.Slice((*reads), func(i, j int) bool { return sorter.Less((*reads)[i], (*reads)[j]) })
}

// QuantizeReads bins quality scores to four Illumina levels. This is lossy —
// the original quality values cannot be recovered. The four output levels are
// Q2, Q11, Q25, Q37 (Phred+33 encoded as '#', ',', ':', 'F'). Reducing the
// quality alphabet from ~40 values to 4 dramatically lowers entropy in the
// quality field, which is often the largest contributor to compressed FASTQ
// file size.
func QuantizeReads(reads []fastq.FastqRead) {
	for i := range reads {
		qual := reads[i].QualityScores()
		binned := make([]byte, len(qual))
		for j, q := range qual {
			binned[j] = quantizeByte(q)
		}
		reads[i].OverrideQual = binned
	}
}

func quantizeByte(q byte) byte {
	// Phred+33 encoding: map raw score to the nearest of four levels.
	score := int(q) - 33
	switch {
	case score < 6:
		return byte(2 + 33)  // Q2  → '#'
	case score < 15:
		return byte(11 + 33) // Q11 → ','
	case score < 27:
		return byte(25 + 33) // Q25 → ':'
	default:
		return byte(37 + 33) // Q37 → 'F'
	}
}
