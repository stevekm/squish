package sort

import (
	"bytes"
	"container/heap"
	fastq "squish/fastq"
)

// FastqReadHeap is a wrapper for FastqRead slice to implement heap.Interface
type FastqReadHeap struct {
	reads          []fastq.FastqRead
	comparisonLess func(a fastq.FastqRead, b fastq.FastqRead) bool
}

func (h FastqReadHeap) Len() int { return len(h.reads) }
func (h FastqReadHeap) Less(i, j int) bool {
	return h.comparisonLess(h.reads[i], h.reads[j])
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
func HeapSort(reads []fastq.FastqRead, comparisonLess func(a fastq.FastqRead, b fastq.FastqRead) bool) {
	heapReads := append([]fastq.FastqRead(nil), reads...)
	h := &FastqReadHeap{reads: heapReads, comparisonLess: comparisonLess}
	heap.Init(h)
	for i := 0; h.Len() > 0; i++ {
		reads[i] = heap.Pop(h).(fastq.FastqRead)
	}
}

func HeapSortSequence(reads *[]fastq.FastqRead) {
	HeapSort(*reads, func(a fastq.FastqRead, b fastq.FastqRead) bool {
		return bytes.Compare(a.Sequence(), b.Sequence()) < 0
	})
}

func HeapSortGC(reads *[]fastq.FastqRead) {
	HeapSort(*reads, func(a fastq.FastqRead, b fastq.FastqRead) bool {
		return a.GCContent < b.GCContent
	})
}
