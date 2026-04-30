package sort

import (
	"bytes"
	"hash/fnv"

	fastq "squish/fastq"
)

// SortStrategy defines the logical ordering for reads.
//
// The in-memory sorter and the external bucket sorter both use this interface,
// which keeps "how reads are compared" independent from "where reads are
// stored while sorting".
type SortStrategy interface {
	Name() string
	Less(a fastq.FastqRead, b fastq.FastqRead) bool
}

// BucketStrategy assigns reads to temporary buckets for external sorting.
//
// OrderedFor reports whether iterating bucket IDs from low to high preserves
// the global ordering for a given sorter. Ordered bucket strategies can be
// sorted one bucket at a time and concatenated. Unordered strategies, such as
// hash buckets, are deterministic and useful for clustering, but do not produce
// a strict global sort for arbitrary sorters without a later merge phase.
type BucketStrategy interface {
	Name() string
	BucketID(read fastq.FastqRead) int
	BucketCount() int
	OrderedFor(sorter SortStrategy) bool
}

// AlphaSort orders reads by sequence bytes.
type AlphaSort struct{}

func (AlphaSort) Name() string { return "alpha" }

func (AlphaSort) Less(a fastq.FastqRead, b fastq.FastqRead) bool {
	if c := bytes.Compare(a.Sequence(), b.Sequence()); c != 0 {
		return c < 0
	}
	// Tie-break on input order so memory and external sorts are reproducible
	// when multiple reads have identical sort keys.
	return a.I < b.I
}

// GCSort orders reads by precomputed GC content.
type GCSort struct{}

func (GCSort) Name() string { return "gc" }

func (GCSort) Less(a fastq.FastqRead, b fastq.FastqRead) bool {
	if a.GCContent == b.GCContent {
		// Keep equal-GC reads deterministic across bucket boundaries.
		return a.I < b.I
	}
	return a.GCContent < b.GCContent
}

// QualitySort orders reads by quality score bytes.
type QualitySort struct{}

func (QualitySort) Name() string { return "qual" }

func (QualitySort) Less(a fastq.FastqRead, b fastq.FastqRead) bool {
	if c := bytes.Compare(a.QualityScores(), b.QualityScores()); c != 0 {
		return c < 0
	}
	// Equal quality strings keep their original input order.
	return a.I < b.I
}

// ClumpSort reuses the clump comparator, which is intended to group similar
// reads for better compression rather than to model a biological ordering.
type ClumpSort struct {
	// K controls the minimizer length used to group reads. A zero value falls
	// back to DefaultClumpKmerLen so older call sites remain valid.
	K int
	// MinCount, when > 1, filters out pivot k-mers that appear fewer than
	// MinCount times across the reads being sorted. Singleton k-mers are
	// usually sequencing errors; excluding them avoids tiny isolated clumps
	// and concentrates reads around confirmed genomic positions.
	MinCount int
	// RComp, when true, reverse-complements reads whose pivot k-mer was chosen
	// from the minus strand. This normalises orientation within each clump so
	// consecutive sequence lines are more byte-similar, increasing LZ77 density.
	RComp bool
	// RawPivot, when true, selects the lex-maximum canonical k-mer instead of
	// the max-hash k-mer. This clusters reads by nucleotide composition (GC
	// content) rather than by hash, which can slightly improve compression for
	// datasets with strong GC bias.
	RawPivot bool
	// Border is the number of bases excluded from each end of the read when
	// selecting the pivot k-mer. Read ends are more error-prone; excluding them
	// avoids pivot k-mers that contain sequencing errors. Clumpify defaults to 1.
	Border int
}

func (ClumpSort) Name() string { return "clump" }

// Less exists only to satisfy SortStrategy. ClumpSort cannot implement a
// correct O(1) comparison without precomputed keys, so all sort paths bypass
// this method: the memory engine calls SortReadsClumpOpts directly, and the
// external engine type-asserts to Sort(). The panic catches any future code
// that accidentally routes clump through a generic Less-based sort loop,
// which would silently recompute the minimizer on every comparison.
func (s ClumpSort) Less(a fastq.FastqRead, b fastq.FastqRead) bool {
	panic("ClumpSort.Less must not be used in a sort loop; call Sort() to get O(n log n) key precomputation")
}

func (s ClumpSort) Sort(reads []fastq.FastqRead) {
	SortReadsClumpOpts(&reads, ClumpSortOptions{
		K:        s.k(),
		MinCount: s.MinCount,
		RComp:    s.RComp,
		RawPivot: s.RawPivot,
		Border:   s.Border,
	})
}

func (s ClumpSort) k() int {
	if s.K < 1 {
		return DefaultClumpKmerLen
	}
	return s.K
}

type BytePrefixBuckets struct {
	name       string
	field      string
	prefixLen  int
	bucketSize int
}

// NewSequencePrefixBuckets creates lexicographically ordered buckets using the
// first prefixLen bytes of the read sequence.
func NewSequencePrefixBuckets(prefixLen int) BytePrefixBuckets {
	return newBytePrefixBuckets("sequence-prefix", "sequence", prefixLen)
}

// NewQualityPrefixBuckets creates lexicographically ordered buckets using the
// first prefixLen bytes of the read quality string.
func NewQualityPrefixBuckets(prefixLen int) BytePrefixBuckets {
	return newBytePrefixBuckets("quality-prefix", "quality", prefixLen)
}

func newBytePrefixBuckets(name string, field string, prefixLen int) BytePrefixBuckets {
	if prefixLen < 1 {
		prefixLen = 1
	}
	bucketSize := 1
	for i := 0; i < prefixLen; i++ {
		bucketSize *= 256
	}
	return BytePrefixBuckets{name: name, field: field, prefixLen: prefixLen, bucketSize: bucketSize}
}

func (b BytePrefixBuckets) Name() string { return b.name }

func (b BytePrefixBuckets) BucketCount() int { return b.bucketSize }

func (b BytePrefixBuckets) BucketID(read fastq.FastqRead) int {
	var key []byte
	switch b.field {
	case "quality":
		key = read.QualityScores()
	default:
		key = read.Sequence()
	}

	bucketID := 0
	for i := 0; i < b.prefixLen; i++ {
		bucketID *= 256
		if i < len(key) {
			// Treat the prefix as a big-endian base-256 integer. This makes
			// bucket ID order match bytewise lexicographic order.
			bucketID += int(key[i])
		}
	}
	return bucketID
}

func (b BytePrefixBuckets) OrderedFor(sorter SortStrategy) bool {
	// Prefix buckets are globally ordered only when the bucketed field matches
	// the sort key. For example, sequence-prefix works for alpha sort but not
	// for quality sort.
	return (b.field == "sequence" && sorter.Name() == "alpha") ||
		(b.field == "quality" && sorter.Name() == "qual")
}

// GCRangeBuckets divides the [0.0, 1.0] GC-content interval into fixed ranges.
type GCRangeBuckets struct {
	bucketCount int
}

// NewGCRangeBuckets creates ordered GC buckets. More buckets reduce the maximum
// bucket size but create more temporary files.
func NewGCRangeBuckets(bucketCount int) GCRangeBuckets {
	if bucketCount < 1 {
		bucketCount = 1
	}
	return GCRangeBuckets{bucketCount: bucketCount}
}

func (b GCRangeBuckets) Name() string { return "gc-range" }

func (b GCRangeBuckets) BucketCount() int { return b.bucketCount }

func (b GCRangeBuckets) BucketID(read fastq.FastqRead) int {
	bucketID := int(read.GCContent * float64(b.bucketCount))
	if bucketID < 0 {
		return 0
	}
	if bucketID >= b.bucketCount {
		return b.bucketCount - 1
	}
	return bucketID
}

func (b GCRangeBuckets) OrderedFor(sorter SortStrategy) bool {
	return sorter.Name() == "gc"
}

// HashBuckets spreads records across a fixed number of buckets using sequence
// and quality bytes. It is useful for bounded memory and clustering, but bucket
// number has no global sort meaning.
type HashBuckets struct {
	name        string
	bucketCount int
	keyFunc     func(fastq.FastqRead) []byte
}

func NewHashBuckets(bucketCount int) HashBuckets {
	return newHashBuckets("hash", bucketCount, func(read fastq.FastqRead) []byte {
		// Use both sequence and quality bytes so hash buckets spread reads more
		// evenly than sequence-only hashing when duplicate reads are common.
		key := make([]byte, 0, len(read.Sequence())+len(read.QualityScores()))
		key = append(key, read.Sequence()...)
		key = append(key, read.QualityScores()...)
		return key
	})
}

func NewClumpBuckets(bucketCount int, k int) HashBuckets {
	// External clump mode hashes the minimizer key rather than the whole record.
	// That puts reads with the same clump key in the same bucket before each
	// bucket is sorted internally.
	return newHashBuckets("clump-minimizer", bucketCount, func(read fastq.FastqRead) []byte {
		return ClumpKeyK(read, k)
	})
}

func newHashBuckets(name string, bucketCount int, keyFunc func(fastq.FastqRead) []byte) HashBuckets {
	if bucketCount < 1 {
		bucketCount = 1
	}
	return HashBuckets{name: name, bucketCount: bucketCount, keyFunc: keyFunc}
}

func (b HashBuckets) Name() string { return b.name }

func (b HashBuckets) BucketCount() int { return b.bucketCount }

func (b HashBuckets) BucketID(read fastq.FastqRead) int {
	h := fnv.New32a()
	h.Write(b.keyFunc(read))
	return int(h.Sum32() % uint32(b.bucketCount))
}

func (b HashBuckets) OrderedFor(sorter SortStrategy) bool {
	// Clump sorting is a compression-oriented grouping pass, so deterministic
	// hash buckets are acceptable even though they are not globally ordered.
	return sorter.Name() == "clump"
}

// StrategyForName maps CLI sort method names to reusable strategy objects.
func StrategyForName(name string) (SortStrategy, bool) {
	switch name {
	case "alpha":
		return AlphaSort{}, true
	case "gc":
		return GCSort{}, true
	case "qual":
		return QualitySort{}, true
	case "clump":
		return ClumpSort{K: DefaultClumpKmerLen}, true
	default:
		return nil, false
	}
}

// DefaultBucketStrategy chooses the most natural external bucket layout for a
// sort strategy. These defaults favor correctness first: alpha, GC, and quality
// get ordered buckets; clump gets deterministic hash buckets.
func DefaultBucketStrategy(sorter SortStrategy, bucketCount int) BucketStrategy {
	switch sorter.Name() {
	case "alpha":
		return NewSequencePrefixBuckets(1)
	case "gc":
		return NewGCRangeBuckets(bucketCount)
	case "qual":
		return NewQualityPrefixBuckets(1)
	case "clump":
		// Preserve the configured clump k-mer length when auto-selecting the
		// external bucket strategy for clump sort.
		clumpSorter, ok := sorter.(ClumpSort)
		if !ok {
			return NewClumpBuckets(bucketCount, DefaultClumpKmerLen)
		}
		return NewClumpBuckets(bucketCount, clumpSorter.k())
	default:
		return NewHashBuckets(bucketCount)
	}
}
