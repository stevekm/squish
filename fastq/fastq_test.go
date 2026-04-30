package fastq

import (
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"

	_io "squish/fastqio"
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

func TestReorderReadsByOrderAppliesPrimaryReadOrder(t *testing.T) {
	dir := t.TempDir()
	orderPath := filepath.Join(dir, "order.txt")
	r2Path := filepath.Join(dir, "r2.fastq")
	outputPath := filepath.Join(dir, "r2.sorted.fastq.gz")

	if err := os.WriteFile(orderPath, []byte("2\n1\n3\n"), 0644); err != nil {
		t.Fatalf("write order: %v", err)
	}
	input := "" +
		"@pair1/2\nTTTT\n+\n!!!!\n" +
		"@pair2/2\nCCCC\n+\n####\n" +
		"@pair3/2\nAAAA\n+\nIIII\n"
	if err := os.WriteFile(r2Path, []byte(input), 0644); err != nil {
		t.Fatalf("write r2: %v", err)
	}

	referenceNames := []string{"pair1", "pair2", "pair3"}
	stats, err := ReorderReadsByOrder(r2Path, outputPath, orderPath, '\n', 3, referenceNames)
	if err != nil {
		t.Fatalf("reorder reads: %v", err)
	}
	if stats.Reads != 3 {
		t.Fatalf("reads = %d, want 3", stats.Reads)
	}

	want := "" +
		"@pair2/2\nCCCC\n+\n####\n" +
		"@pair1/2\nTTTT\n+\n!!!!\n" +
		"@pair3/2\nAAAA\n+\nIIII\n"
	if got := readGzipFile(t, outputPath); got != want {
		t.Fatalf("reordered output = %q, want %q", got, want)
	}
}

func TestReorderReadsByOrderRejectsReadCountMismatch(t *testing.T) {
	dir := t.TempDir()
	orderPath := filepath.Join(dir, "order.txt")
	r2Path := filepath.Join(dir, "r2.fastq")
	outputPath := filepath.Join(dir, "r2.sorted.fastq.gz")

	if err := os.WriteFile(orderPath, []byte("1\n2\n"), 0644); err != nil {
		t.Fatalf("write order: %v", err)
	}
	if err := os.WriteFile(r2Path, []byte("@pair1/2\nTTTT\n+\n!!!!\n"), 0644); err != nil {
		t.Fatalf("write r2: %v", err)
	}

	if _, err := ReorderReadsByOrder(r2Path, outputPath, orderPath, '\n', 2, nil); err == nil {
		t.Fatalf("expected read count mismatch error")
	}
}

func TestReorderReadsByOrderRejectsReadNameMismatch(t *testing.T) {
	dir := t.TempDir()
	orderPath := filepath.Join(dir, "order.txt")
	r2Path := filepath.Join(dir, "r2.fastq")
	outputPath := filepath.Join(dir, "r2.sorted.fastq.gz")

	if err := os.WriteFile(orderPath, []byte("1\n"), 0644); err != nil {
		t.Fatalf("write order: %v", err)
	}
	if err := os.WriteFile(r2Path, []byte("@wrong/2\nTTTT\n+\n!!!!\n"), 0644); err != nil {
		t.Fatalf("write r2: %v", err)
	}

	if _, err := ReorderReadsByOrder(r2Path, outputPath, orderPath, '\n', 1, []string{"pair1"}); err == nil {
		t.Fatalf("expected read name mismatch error")
	}
}

func readGzipFile(t *testing.T, path string) string {
	t.Helper()

	file, err := os.Open(path)
	if err != nil {
		t.Fatalf("open gzip: %v", err)
	}
	defer file.Close()

	reader, err := gzip.NewReader(file)
	if err != nil {
		t.Fatalf("create gzip reader: %v", err)
	}
	defer reader.Close()

	data, err := io.ReadAll(reader)
	if err != nil {
		t.Fatalf("read gzip: %v", err)
	}
	return string(data)
}
