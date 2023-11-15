package sort

import (
	go_sort "sort"
	fastq "squish/fastq"
)

func SortReadsSequence(reads *[]fastq.FastqRead) {
	// put the read sorting logic in here!
	//
	// sort in-place based on the string of the sequence
	go_sort.Slice((*reads), func(i, j int) bool { return string((*reads)[i].Sequence) < string((*reads)[j].Sequence) })
}

func SortReadsGC(reads *[]fastq.FastqRead) {
	go_sort.Slice((*reads), func(i, j int) bool { return (*reads)[i].GCContent < (*reads)[j].GCContent })
}

func SortReadsQual(reads *[]fastq.FastqRead) {
	go_sort.Slice((*reads), func(i, j int) bool { return string((*reads)[i].QualityScores) < string((*reads)[j].QualityScores) })
}
