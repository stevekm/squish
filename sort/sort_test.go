package sort

import (
	"os"
	"path/filepath"
	fastq "squish/fastq"
	_io "squish/fastqio"
	"testing"
)

func TestSortReadsSequenceExample(t *testing.T) {
	input := "" +
		"@read_t\nTTTT\n+\nIIII\n" +
		"@read_a\nAAAA\n+\n####\n" +
		"@read_c\nCCCC\n+\n!!!!\n"
	reads := loadReadsFromString(t, input)

	SortReadsSequence(&reads)

	want := []string{
		"@read_a\nAAAA\n+\n####\n",
		"@read_c\nCCCC\n+\n!!!!\n",
		"@read_t\nTTTT\n+\nIIII\n",
	}
	assertRecords(t, reads, want)
	assertExternalSortOutput(t, input, AlphaSort{}, NewSequencePrefixBuckets(1), want)
}

func TestSortReadsGCExample(t *testing.T) {
	input := "" +
		"@all_gc\nGGGG\n+\nIIII\n" +
		"@no_gc\nATAT\n+\nIIII\n" +
		"@half_gc\nACAT\n+\nIIII\n"
	reads := loadReadsFromString(t, input)

	SortReadsGC(&reads)

	want := []string{
		"@no_gc\nATAT\n+\nIIII\n",
		"@half_gc\nACAT\n+\nIIII\n",
		"@all_gc\nGGGG\n+\nIIII\n",
	}
	assertRecords(t, reads, want)
	assertExternalSortOutput(t, input, GCSort{}, NewGCRangeBuckets(4), want)
}

func TestSortReadsQualExample(t *testing.T) {
	input := "" +
		"@high_quality\nAAAA\n+\nIIII\n" +
		"@low_quality\nCCCC\n+\n!!!!\n" +
		"@mid_quality\nGGGG\n+\n####\n"
	reads := loadReadsFromString(t, input)

	SortReadsQual(&reads)

	want := []string{
		"@low_quality\nCCCC\n+\n!!!!\n",
		"@mid_quality\nGGGG\n+\n####\n",
		"@high_quality\nAAAA\n+\nIIII\n",
	}
	assertRecords(t, reads, want)
	assertExternalSortOutput(t, input, QualitySort{}, NewQualityPrefixBuckets(1), want)
}

func TestSortReadsClumpExample(t *testing.T) {
	input := "" +
		"@poly_g\nGGGGGG\n+\nIIIIII\n" +
		"@has_aaa\nTTTAAA\n+\nIIIIII\n" +
		"@poly_c\nCCCCCC\n+\nIIIIII\n" +
		"@starts_aaa\nAAAGGG\n+\nIIIIII\n"
	reads := loadReadsFromString(t, input)

	SortReadsClumpK(&reads, 3)

	want := []string{
		"@starts_aaa\nAAAGGG\n+\nIIIIII\n",
		"@has_aaa\nTTTAAA\n+\nIIIIII\n",
		"@poly_c\nCCCCCC\n+\nIIIIII\n",
		"@poly_g\nGGGGGG\n+\nIIIIII\n",
	}
	assertRecords(t, reads, want)
	assertExternalSortOutput(t, input, ClumpSort{K: 3}, NewClumpBuckets(1, 3), want)
}

func loadReadsFromString(t *testing.T, input string) []fastq.FastqRead {
	t.Helper()

	inputPath := filepath.Join(t.TempDir(), "input.fastq")
	if err := os.WriteFile(inputPath, []byte(input), 0644); err != nil {
		t.Fatalf("write input fastq: %v", err)
	}

	reader := _io.GetReader(inputPath)
	defer reader.Close()

	delim := byte('\n')
	reads := []fastq.FastqRead{}
	fastq.LoadReads(&reads, reader, &delim)
	return reads
}

func assertRecords(t *testing.T, reads []fastq.FastqRead, want []string) {
	t.Helper()

	if len(reads) != len(want) {
		t.Fatalf("read count = %d, want %d", len(reads), len(want))
	}
	for i, read := range reads {
		if got := string(read.Record()); got != want[i] {
			t.Fatalf("record %d = %q, want %q", i, got, want[i])
		}
	}
}

func assertExternalSortOutput(t *testing.T, input string, sorter SortStrategy, bucketer BucketStrategy, want []string) {
	t.Helper()

	dir := t.TempDir()
	inputPath := filepath.Join(dir, "input.fastq")
	outputPath := filepath.Join(dir, "output.fastq.gz")
	orderPath := filepath.Join(dir, "order.txt")
	tempDir := filepath.Join(dir, "tmp")
	if err := os.WriteFile(inputPath, []byte(input), 0644); err != nil {
		t.Fatalf("write external input fastq: %v", err)
	}

	config := ExternalBucketConfig{
		InputFilepath:  inputPath,
		OutputFilepath: outputPath,
		OrderFilepath:  orderPath,
		TempDir:        tempDir,
		RecordDelim:    '\n',
	}
	stats, err := RunExternalBucketSort(config, sorter, bucketer)
	if err != nil {
		t.Fatalf("external %s sort with %s buckets: %v", sorter.Name(), bucketer.Name(), err)
	}
	if stats.Reads != len(want) {
		t.Fatalf("external reads = %d, want %d", stats.Reads, len(want))
	}

	got := readGzipFile(t, outputPath)
	wantOutput := joinRecords(want)
	if got != wantOutput {
		t.Fatalf("external output = %q, want %q", got, wantOutput)
	}
}

func joinRecords(records []string) string {
	output := ""
	for _, record := range records {
		output += record
	}
	return output
}
