package sort

import (
	"bytes"
	go_sort "sort"
	fastq "squish/fastq"
)

// using built-in Go sort methods
func SortReadsSequence(reads *[]fastq.FastqRead) {
	// put the read sorting logic in here!
	//
	// sort in-place based on the string of the sequence
	go_sort.Slice((*reads), func(i, j int) bool {
		return bytes.Compare((*reads)[i].Sequence(), (*reads)[j].Sequence()) < 0
	})
}

func SortReadsGC(reads *[]fastq.FastqRead) {
	go_sort.Slice((*reads), func(i, j int) bool { return (*reads)[i].GCContent < (*reads)[j].GCContent })
}

func SortReadsQual(reads *[]fastq.FastqRead) {
	go_sort.Slice((*reads), func(i, j int) bool {
		return bytes.Compare((*reads)[i].QualityScores(), (*reads)[j].QualityScores()) < 0
	})
}
