package sort

import (
	"testing"

	fastq "squish/fastq"
)

// makeClumpRead builds a minimal clumpRead for unit-testing ClumpReadLess.
// The embedded FastqRead has a real arena so Sequence(), QualityScores(),
// and Id() all work without panicking.
func makeClumpRead(seq, qual, id string, i int, key []byte, pivotPos int) clumpRead {
	header := "@" + id + "\n"
	plus := "+\n"
	data := header + seq + "\n" + plus + qual + "\n"
	arena := &fastq.FastqArena{Data: []byte(data)}
	read := fastq.FastqRead{
		Arena:              arena,
		RecordOffset:       0,
		RecordSize:         len(data),
		IdOffset:           0,
		IdSize:             len(header),
		SequenceOffset:     len(header),
		SequenceSize:       len(seq) + 1,
		PlusOffset:         len(header) + len(seq) + 1,
		PlusSize:           len(plus),
		QualityScoreOffset: len(header) + len(seq) + 1 + len(plus),
		QualityScoreSize:   len(qual) + 1,
		I:                  i,
	}
	return clumpRead{read: read, key: key, pivotPos: pivotPos}
}

// reverseComplement helpers

func TestReverseComplementAllBases(t *testing.T) {
	cases := []struct{ in, want string }{
		{"A", "T"},
		{"C", "G"},
		{"G", "C"},
		{"T", "A"},
		{"N", "N"},
		{"ACGT", "ACGT"}, // palindrome
		{"AAAA", "TTTT"},
		{"TTTT", "AAAA"},
		{"CCCC", "GGGG"},
		{"GCGC", "GCGC"}, // palindrome
	}
	for _, tc := range cases {
		got := string(reverseComplement([]byte(tc.in)))
		if got != tc.want {
			t.Errorf("reverseComplement(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

// clumpMinimizerFull

func TestClumpMinimizerFullRCFlipped(t *testing.T) {
	// "TTTTTT" k=2: every k-mer is TT whose canonical form is AA (A < T).
	// The pivot is on the minus strand so rcFlipped should be true.
	key, _, rcFlipped := clumpMinimizerFull([]byte("TTTTTT"), 2, false, nil, 0)
	if string(key) != "AA" {
		t.Fatalf("key = %q, want AA", key)
	}
	if !rcFlipped {
		t.Fatal("rcFlipped should be true for all-T sequence")
	}
}

func TestClumpMinimizerFullForwardStrand(t *testing.T) {
	// "AAAAAA" k=2: every k-mer is AA, which is already canonical (A < T).
	// rcFlipped should be false.
	key, _, rcFlipped := clumpMinimizerFull([]byte("AAAAAA"), 2, false, nil, 0)
	if string(key) != "AA" {
		t.Fatalf("key = %q, want AA", key)
	}
	if rcFlipped {
		t.Fatal("rcFlipped should be false for all-A sequence")
	}
}

func TestClumpMinimizerFullBorderExcludesEnds(t *testing.T) {
	// "TTTAAA" k=3 border=1:
	//   valid positions 1 and 2 only (border=1 excludes pos 0 and pos 3)
	//   pos 1: "TTA" → RC "TAA" (A < T so TAA is canonical), rcFlipped=true
	//   pos 2: "TAA" → RC "TTA" (TAA < TTA so TAA is canonical), rcFlipped=false
	//   both give canonical TAA, so the key must be "TAA".
	key, _, _ := clumpMinimizerFull([]byte("TTTAAA"), 3, false, nil, 1)
	if string(key) != "TAA" {
		t.Fatalf("border=1 key = %q, want TAA", key)
	}
}

func TestClumpMinimizerFullShortSeq(t *testing.T) {
	// Sequence "ACG" is shorter than k=4: the whole sequence is the single k-mer.
	// RC of "ACG" is "CGT"; "ACG" < "CGT" so forward is canonical.
	key, pos, rcFlipped := clumpMinimizerFull([]byte("ACG"), 4, false, nil, 0)
	if string(key) != "ACG" {
		t.Fatalf("key = %q, want ACG", key)
	}
	if pos != 0 {
		t.Fatalf("pos = %d, want 0", pos)
	}
	if rcFlipped {
		t.Fatal("rcFlipped should be false (ACG < CGT)")
	}
}

func TestClumpMinimizerFullNoEligible(t *testing.T) {
	// When every k-mer is rejected by the eligible filter, nil key is returned.
	key, _, _ := clumpMinimizerFull([]byte("ACGTACGT"), 3, false, func([]byte) bool { return false }, 0)
	if key != nil {
		t.Fatalf("expected nil key when no k-mers eligible, got %q", key)
	}
}

// ClumpReadLess

func TestClumpReadLessByKey(t *testing.T) {
	a := makeClumpRead("AAAA", "IIII", "a", 1, []byte("AAA"), 0)
	b := makeClumpRead("TTTT", "IIII", "b", 2, []byte("TTT"), 0)
	if !ClumpReadLess(a, b) {
		t.Fatal("expected key AAA < TTT")
	}
	if ClumpReadLess(b, a) {
		t.Fatal("expected key TTT not < AAA")
	}
}

func TestClumpReadLessNilKeyBeforeNonNil(t *testing.T) {
	// Nil key (unclustered reads) sorts before any real key.
	a := makeClumpRead("AAAA", "IIII", "a", 1, nil, 0)
	b := makeClumpRead("TTTT", "IIII", "b", 2, []byte("TTT"), 0)
	if !ClumpReadLess(a, b) {
		t.Fatal("expected nil key < non-nil key")
	}
}

func TestClumpReadLessByPivotPos(t *testing.T) {
	// Same key, different pivot positions: smaller position sorts first.
	key := []byte("ACG")
	a := makeClumpRead("ACGAAA", "IIIIII", "a", 1, key, 0)
	b := makeClumpRead("AAAACG", "IIIIII", "b", 2, key, 3)
	if !ClumpReadLess(a, b) {
		t.Fatal("expected pivotPos=0 < pivotPos=3")
	}
	if ClumpReadLess(b, a) {
		t.Fatal("expected pivotPos=3 not < pivotPos=0")
	}
}

