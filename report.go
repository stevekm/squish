package squish

import (
	"encoding/json"
	"fmt"
	"log/slog"
	"os"
	"strings"
)

type FileReport struct {
	Path       string `json:"path"`
	Argument   string `json:"argument,omitempty"`
	SizeBytes  int64  `json:"size_bytes,omitempty"`
	SizeHuman  string `json:"size_human,omitempty"`
	Descriptor string `json:"descriptor,omitempty"`
}

type ProfileReport struct {
	Directory string `json:"directory"`
	CPUPath   string `json:"cpu_path"`
	MemPath   string `json:"mem_path"`
}

type BucketReport struct {
	Strategy   string `json:"strategy"`
	Count      int    `json:"count"`
	Used       int    `json:"used"`
	TempDir    string `json:"temp_dir,omitempty"`
	OrderedFor bool   `json:"ordered_for_sorter"`
}

type PairedReport struct {
	Input             FileReport `json:"input"`
	Output            FileReport `json:"output"`
	Reads             int        `json:"reads"`
	UncompressedBytes int        `json:"uncompressed_bytes"`
	OutputSizeBytes   int64      `json:"output_size_bytes"`
}

type Report struct {
	Version              string         `json:"version"`
	StartedAt            string         `json:"started_at"`
	FinishedAt           string         `json:"finished_at"`
	Duration             string         `json:"duration"`
	DurationMilliseconds int64          `json:"duration_ms"`
	SortMethod           string         `json:"sort_method"`
	SortDescription      string         `json:"sort_description"`
	SortEngine           string         `json:"sort_engine"`
	ClumpKmerLength      int            `json:"clump_kmer_length"`
	Input                FileReport     `json:"input"`
	Output               FileReport     `json:"output"`
	OrderFile            FileReport     `json:"order_file"`
	ReportFile           FileReport     `json:"report_file"`
	ManifestFile         FileReport     `json:"manifest_file"`
	PairedOutputs        []PairedReport `json:"paired_outputs,omitempty"`
	Profile              ProfileReport  `json:"profile"`
	Bucket               *BucketReport  `json:"bucket,omitempty"`
	Reads                int            `json:"reads"`
	UncompressedBytes    int            `json:"uncompressed_bytes"`
	OutputSizeBytes      int64          `json:"output_size_bytes"`
	SizeDifferenceBytes  int64          `json:"size_difference_bytes"`
	CompressionRatio     float64        `json:"compression_ratio"`
	SizeReductionRatio   float64        `json:"size_reduction_ratio"`
}

func WriteReport(report Report, reportPath string) error {
	reportJSON, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		return fmt.Errorf("marshal report: %w", err)
	}
	if err := os.WriteFile(reportPath, reportJSON, 0644); err != nil {
		return fmt.Errorf("write report file %q: %w", reportPath, err)
	}
	slog.Debug("report written", "path", reportPath)
	return nil
}

func WriteManifest(outputPaths []string, manifestPath string) error {
	lines := make([]string, 0, len(outputPaths))
	for _, outputPath := range outputPaths {
		absolutePath, err := AbsolutePath(outputPath)
		if err != nil {
			return err
		}
		lines = append(lines, absolutePath)
	}
	manifestText := strings.Join(lines, "\n")
	if manifestText != "" {
		manifestText += "\n"
	}
	if err := os.WriteFile(manifestPath, []byte(manifestText), 0644); err != nil {
		return fmt.Errorf("write manifest file %q: %w", manifestPath, err)
	}
	slog.Debug("manifest written", "path", manifestPath, "outputs", len(lines))
	return nil
}
