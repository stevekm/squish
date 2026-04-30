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
