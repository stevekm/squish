package fastq

import (
	"os"
	"path/filepath"
	"testing"

	_io "squish/io"
)

func TestLoadReadsStoresRecordsInArena(t *testing.T) {
	input := []byte("@read2\nTTTT\n+\n!!!!\n@read1\nACGT\n+\nIIII\n")
	inputPath := filepath.Join(t.TempDir(), "reads.fastq")
	if err := os.WriteFile(inputPath, input, 0644); err != nil {
		t.Fatalf("write input fastq: %v", err)
	}

	reader := _io.GetReader(inputPath)
	defer reader.Close()

	delim := byte('\n')
	reads := []FastqRead{}
	totalSize := LoadReads(&reads, reader, &delim)

	if totalSize != len(input) {
		t.Fatalf("total size = %d, want %d", totalSize, len(input))
	}
	if len(reads) != 2 {
		t.Fatalf("read count = %d, want 2", len(reads))
	}
	if reads[0].Arena == nil || reads[0].Arena != reads[1].Arena {
		t.Fatalf("reads should share the same arena")
	}
	if got := string(reads[0].Record()); got != "@read2\nTTTT\n+\n!!!!\n" {
		t.Fatalf("first record = %q", got)
	}
	if got := string(reads[1].Sequence()); got != "ACGT\n" {
		t.Fatalf("second sequence = %q", got)
	}
	if got := string(reads[1].QualityScores()); got != "IIII\n" {
		t.Fatalf("second quality = %q", got)
	}
}
