package sort

import (
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"
)

func TestRunExternalBucketSortAlpha(t *testing.T) {
	dir := t.TempDir()
	inputPath := filepath.Join(dir, "input.fastq")
	outputPath := filepath.Join(dir, "output.fastq.gz")
	orderPath := filepath.Join(dir, "order.txt")
	tempDir := filepath.Join(dir, "tmp")

	input := []byte("@read2\nTTTT\n+\n!!!!\n@read1\nACGT\n+\nIIII\n")
	if err := os.WriteFile(inputPath, input, 0644); err != nil {
		t.Fatalf("write input: %v", err)
	}

	config := ExternalBucketConfig{
		InputFilepath:  inputPath,
		OutputFilepath: outputPath,
		OrderFilepath:  orderPath,
		TempDir:        tempDir,
		RecordDelim:    '\n',
	}
	if err := RunExternalBucketSort(config, AlphaSort{}, NewSequencePrefixBuckets(1)); err != nil {
		t.Fatalf("external sort: %v", err)
	}

	got := readGzipFile(t, outputPath)
	want := "@read1\nACGT\n+\nIIII\n@read2\nTTTT\n+\n!!!!\n"
	if got != want {
		t.Fatalf("output = %q, want %q", got, want)
	}

	order, err := os.ReadFile(orderPath)
	if err != nil {
		t.Fatalf("read order: %v", err)
	}
	if got, want := string(order), "2\n1\n"; got != want {
		t.Fatalf("order = %q, want %q", got, want)
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

	data, err := os.ReadFile(path)
	if err == nil && len(data) == 0 {
		t.Fatalf("gzip file is empty")
	}

	uncompressed, err := io.ReadAll(reader)
	if err != nil {
		t.Fatalf("read gzip: %v", err)
	}
	return string(uncompressed)
}
