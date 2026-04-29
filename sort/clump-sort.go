package sort

import (
	"bytes"
	go_sort "sort"
	fastq "squish/fastq"
)

func SortReadsClump(reads *[]fastq.FastqRead) {
	const prefixLen = 16

	// Group reads by a simple clump key based on the first few bases of the sequence
	// and the first few quality characters. This creates buckets of similar reads.
	buckets := make(map[string][]fastq.FastqRead, len(*reads)/8)
	bucketKeys := make([]string, 0, len(*reads))

	for _, read := range *reads {
		key := clumpBucketKey(&read, prefixLen)
		if _, ok := buckets[key]; !ok {
			// remember the order of bucket keys so we can process buckets deterministically
			bucketKeys = append(bucketKeys, key)
		}
		buckets[key] = append(buckets[key], read)
	}

	// Sort bucket keys to ensure the bucket ordering is stable across runs.
	go_sort.Strings(bucketKeys)

	// Sort each bucket internally by the stronger full-read comparison.
	sortedReads := make([]fastq.FastqRead, 0, len(*reads))
	for _, key := range bucketKeys {
		bucket := buckets[key]
		go_sort.Slice(bucket, func(i, j int) bool {
			return ClumpCompare(bucket[i], bucket[j])
		})
		sortedReads = append(sortedReads, bucket...)
	}

	// Replace the original read slice with the new clump-ordered slice.
	*reads = sortedReads
}

func clumpBucketKey(read *fastq.FastqRead, prefixLen int) string {
	seq := read.Sequence()
	if len(seq) > prefixLen {
		seq = seq[:prefixLen]
	}

	qual := read.QualityScores()
	if len(qual) > prefixLen {
		qual = qual[:prefixLen]
	}

	// Use a simple composite key to cluster reads with similar sequence and quality prefixes.
	return string(seq) + "|" + string(qual)
}

func ClumpCompare(a, b fastq.FastqRead) bool {
	// Prefer reads with identical full sequence content.
	if c := bytes.Compare(a.Sequence(), b.Sequence()); c != 0 {
		return c < 0
	}
	// Tie break on quality string similarity.
	if c := bytes.Compare(a.QualityScores(), b.QualityScores()); c != 0 {
		return c < 0
	}
	// Keep headers stable so repeated Id prefixes stay together.
	if c := bytes.Compare(a.Id(), b.Id()); c != 0 {
		return c < 0
	}
	// Finally preserve original file order for identical reads.
	return a.I < b.I
}
