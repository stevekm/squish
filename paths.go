package squish

import (
	"fmt"
	"log/slog"
	"os"
	"path/filepath"
	"strings"

	"code.cloudfoundry.org/bytefmt"
)

func OutputPath(outputDir string, itemPath string) (string, error) {
	outputPath := filepath.Join(outputDir, itemPath)
	relPath, err := filepath.Rel(outputDir, outputPath)
	if err != nil || relPath == ".." || strings.HasPrefix(relPath, ".."+string(os.PathSeparator)) {
		return "", fmt.Errorf("output item path %q must stay within output directory %q", itemPath, outputDir)
	}
	return outputPath, nil
}

func AbsolutePath(path string) (string, error) {
	absolutePath, err := filepath.Abs(path)
	if err != nil {
		return "", fmt.Errorf("resolve absolute path for %q: %w", path, err)
	}
	return absolutePath, nil
}

func SplitPathList(value string) []string {
	items := []string{}
	for _, item := range strings.FieldsFunc(value, func(r rune) bool {
		return r == ',' || r == ';'
	}) {
		item = strings.TrimSpace(item)
		if item != "" {
			items = append(items, item)
		}
	}
	return items
}

func StripFastqExtensions(path string) string {
	base := filepath.Base(path)
	for _, suffix := range []string{".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"} {
		if strings.HasSuffix(base, suffix) {
			return strings.TrimSuffix(base, suffix)
		}
	}
	return base
}

func DerivePairedOutputArg(inputPath string) string {
	return StripFastqExtensions(inputPath) + ".sorted.fastq.gz"
}

func GetFileSize(path string) (int64, error) {
	fi, err := os.Stat(path)
	if err != nil {
		return 0, err
	}
	return fi.Size(), nil
}

func LogFileSize(path string, filetype string) int64 {
	inputFileSize, err := GetFileSize(path)
	if err != nil {
		slog.Debug("could not get size for file", "path", path, "error", err)
	}
	slog.Debug("file size", "type", filetype, "path", path, "size", bytefmt.ByteSize(uint64(inputFileSize)))
	return inputFileSize
}
