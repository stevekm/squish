package sort

import (
	"container/heap"
	go_sort "sort"
	fastq "squish/fastq"
)

// using built-in Go sort methods
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

//
// using custom sort methods
//


// FastqReadHeap is a wrapper for FastqRead slice to implement heap.Interface
type FastqReadHeap struct {
	reads         []fastq.FastqRead
	comparisonKey func(fr fastq.FastqRead) interface{}
}

func (h FastqReadHeap) Len() int { return len(h.reads) }
func (h FastqReadHeap) Less(i, j int) bool {
	// return h.comparisonKey(h.reads[i]) < h.comparisonKey(h.reads[j])
	keyI := h.comparisonKey(h.reads[i])
	keyJ := h.comparisonKey(h.reads[j])

	switch keyI := keyI.(type) {
	case float64:
		return keyI < keyJ.(float64)
	case string:
		return keyI < keyJ.(string)
	default:
		// Handle other types if necessary
		return false
	}
}
func (h FastqReadHeap) Swap(i, j int) { h.reads[i], h.reads[j] = h.reads[j], h.reads[i] }

// Push appends an element to the heap.
func (h *FastqReadHeap) Push(x interface{}) {
	h.reads = append(h.reads, x.(fastq.FastqRead))
}

// Pop removes and returns the minimum element (according to Less) from the heap.
func (h *FastqReadHeap) Pop() interface{} {
	old := h.reads
	n := len(old)
	x := old[n-1]
	h.reads = old[0 : n-1]
	return x
}

// HeapSort sorts the FastqRead slice in-place using heap sort.
func HeapSort(reads []fastq.FastqRead, comparisonKey func(fr fastq.FastqRead) interface{}) {
	h := &FastqReadHeap{reads: reads, comparisonKey: comparisonKey}
	heap.Init(h)
	for h.Len() > 0 {
		heap.Pop(h)
	}
}

func HeapSortSequence(reads *[]fastq.FastqRead) {
	HeapSort(*reads, func(fr fastq.FastqRead) interface{} {
		return string(fr.Sequence)
	})
}

func HeapSortGC(reads *[]fastq.FastqRead) {
	HeapSort(*reads, func(fr fastq.FastqRead) interface{} {
		return fr.GCContent
	})
}

