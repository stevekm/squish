package squish

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestRunExternalAlphaWritesReportManifestAndPairedOutput(t *testing.T) {
	dir := t.TempDir()
	r1Path := filepath.Join(dir, "r1.fastq")
	r2Path := filepath.Join(dir, "r2.fastq")
	outDir := filepath.Join(dir, "out")

	r1 := "" +
		"@pair2/1\nTTTT\n+\n####\n" +
		"@pair1/1\nAAAA\n+\nIIII\n" +
		"@pair3/1\nCCCC\n+\n!!!!\n"
	r2 := "" +
		"@pair2/2\nGGGG\n+\n####\n" +
		"@pair1/2\nTTTT\n+\nIIII\n" +
		"@pair3/2\nAAAA\n+\n!!!!\n"
	if err := os.WriteFile(r1Path, []byte(r1), 0644); err != nil {
		t.Fatalf("write r1: %v", err)
	}
	if err := os.WriteFile(r2Path, []byte(r2), 0644); err != nil {
		t.Fatalf("write r2: %v", err)
	}

	result, err := Run(context.Background(), Config{
		SortMethod:            "alpha",
		SortEngine:            "external",
		BucketStrategy:        "auto",
		BucketCount:           8,
		InputFilepath:         r1Path,
		OutputFilenameArg:     "r1.sorted.fastq.gz",
		OutputDir:             outDir,
		PairedInputArgs:       []string{r2Path},
		PairedInputFilepaths:  []string{r2Path},
		PairedOutputArgs:      []string{"r2.sorted.fastq.gz"},
		PairedOutputFilepaths: []string{filepath.Join(outDir, "r2.sorted.fastq.gz")},
		CheckPairs:            true,
		RecordDelim:           RecordDelim,
		RecordHeaderChar:      FastqHeaderChar,
		OrderFilename:         filepath.Join(outDir, "order.txt"),
		ReportFilename:        filepath.Join(outDir, "report.json"),
		ManifestFilename:      filepath.Join(outDir, "manifest.txt"),
		TempDir:               filepath.Join(outDir, "tmp"),
		ProfileDir:            filepath.Join(outDir, "profile.alpha"),
		CPUProfilePath:        filepath.Join(outDir, "profile.alpha", "cpu.prof"),
		MemProfilePath:        filepath.Join(outDir, "profile.alpha", "mem.prof"),
	})
	if err != nil {
		t.Fatalf("run squish: %v", err)
	}

	if result.Report.Reads != 3 {
		t.Fatalf("report reads = %d, want 3", result.Report.Reads)
	}
	if len(result.Report.PairedOutputs) != 1 {
		t.Fatalf("paired outputs = %d, want 1", len(result.Report.PairedOutputs))
	}

	for _, path := range []string{
		filepath.Join(outDir, "r1.sorted.fastq.gz"),
		filepath.Join(outDir, "r2.sorted.fastq.gz"),
		filepath.Join(outDir, "order.txt"),
		filepath.Join(outDir, "report.json"),
		filepath.Join(outDir, "manifest.txt"),
	} {
		if _, err := os.Stat(path); err != nil {
			t.Fatalf("expected output %s: %v", path, err)
		}
	}

	manifest, err := os.ReadFile(filepath.Join(outDir, "manifest.txt"))
	if err != nil {
		t.Fatalf("read manifest: %v", err)
	}
	manifestText := string(manifest)
	if !strings.Contains(manifestText, filepath.Join(outDir, "r1.sorted.fastq.gz")) {
		t.Fatalf("manifest missing primary output: %q", manifestText)
	}
	if !strings.Contains(manifestText, filepath.Join(outDir, "r2.sorted.fastq.gz")) {
		t.Fatalf("manifest missing paired output: %q", manifestText)
	}
}
