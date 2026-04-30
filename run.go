package squish

import (
	"context"
	"fmt"
	"log/slog"
	"time"

	"code.cloudfoundry.org/bytefmt"
)

func Run(ctx context.Context, config Config) (Result, error) {
	select {
	case <-ctx.Done():
		return Result{}, ctx.Err()
	default:
	}

	config, sortDefinition, err := normalizeConfig(config)
	if err != nil {
		return Result{}, err
	}
	if err := ensureOutputDirs(config); err != nil {
		return Result{}, err
	}

	config.InputFileSize = LogFileSize(config.InputFilepath, "Input")

	cpuFile, memFile, err := startProfiling(config.CPUProfilePath, config.MemProfilePath)
	if err != nil {
		return Result{}, err
	}
	defer func() {
		if err := stopProfiling(cpuFile, memFile); err != nil {
			slog.Debug("could not stop profiling cleanly", "error", err)
		}
	}()

	var runStats RunStats
	switch config.SortEngine {
	case "memory":
		runStats, err = RunSort(config, sortDefinition)
	case "external":
		runStats, err = RunExternalSort(config, sortDefinition)
	default:
		return Result{}, fmt.Errorf("unknown sort engine: %s", config.SortEngine)
	}
	if err != nil {
		return Result{}, err
	}

	pairedStats, err := RunPairedReorders(config, runStats.Reads)
	if err != nil {
		return Result{}, err
	}

	timeStop := time.Now()
	timeDuration := timeStop.Sub(config.TimeStart)
	outputFileSize := LogFileSize(config.OutputFilepath, "Output")
	sizeDifference := config.InputFileSize - outputFileSize
	compressionRatio := 0.0
	sizeReductionRatio := 0.0
	if config.InputFileSize > 0 {
		compressionRatio = float64(outputFileSize) / float64(config.InputFileSize)
		sizeReductionRatio = float64(sizeDifference) / float64(config.InputFileSize)
	}

	var bucketReport *BucketReport
	if config.SortEngine == "external" {
		bucketer, err := GetBucketStrategy(config, sortDefinition)
		if err != nil {
			return Result{}, err
		}
		bucketReport = &BucketReport{
			Strategy:   runStats.BucketName,
			Count:      runStats.BucketCount,
			Used:       runStats.BucketsUsed,
			TempDir:    runStats.BucketTempDir,
			OrderedFor: bucketer.OrderedFor(sortDefinition.Strategy),
		}
	}

	outputAbsolutePath, err := AbsolutePath(config.OutputFilepath)
	if err != nil {
		return Result{}, err
	}
	pairedReports := make([]PairedReport, 0, len(pairedStats))
	manifestOutputPaths := []string{config.OutputFilepath}
	for _, pairedStat := range pairedStats {
		pairedOutputAbsolutePath, err := AbsolutePath(pairedStat.OutputPath)
		if err != nil {
			return Result{}, err
		}
		manifestOutputPaths = append(manifestOutputPaths, pairedStat.OutputPath)
		pairedReports = append(pairedReports, PairedReport{
			Input: FileReport{
				Path:      pairedStat.InputPath,
				Argument:  pairedStat.InputArgument,
				SizeBytes: pairedStat.InputSizeBytes,
				SizeHuman: bytefmt.ByteSize(uint64(pairedStat.InputSizeBytes)),
			},
			Output: FileReport{
				Path:      pairedOutputAbsolutePath,
				Argument:  pairedStat.OutputArgument,
				SizeBytes: pairedStat.OutputSizeBytes,
				SizeHuman: bytefmt.ByteSize(uint64(pairedStat.OutputSizeBytes)),
			},
			Reads:             pairedStat.Reads,
			UncompressedBytes: pairedStat.UncompressedBytes,
			OutputSizeBytes:   pairedStat.OutputSizeBytes,
		})
	}
	if err := WriteManifest(manifestOutputPaths, config.ManifestFilename); err != nil {
		return Result{}, err
	}
	manifestFileSize := LogFileSize(config.ManifestFilename, "Manifest")

	report := Report{
		Version:              Version,
		StartedAt:            config.TimeStart.Format(time.RFC3339Nano),
		FinishedAt:           timeStop.Format(time.RFC3339Nano),
		Duration:             timeDuration.String(),
		DurationMilliseconds: timeDuration.Milliseconds(),
		SortMethod:           sortDefinition.CLIArg,
		SortDescription:      sortDefinition.Description,
		SortEngine:           config.SortEngine,
		ClumpKmerLength:      config.ClumpKmerLen,
		Input: FileReport{
			Path:      config.InputFilepath,
			SizeBytes: config.InputFileSize,
			SizeHuman: bytefmt.ByteSize(uint64(config.InputFileSize)),
		},
		Output: FileReport{
			Path:      outputAbsolutePath,
			Argument:  config.OutputFilenameArg,
			SizeBytes: outputFileSize,
			SizeHuman: bytefmt.ByteSize(uint64(outputFileSize)),
		},
		OrderFile: FileReport{
			Path: config.OrderFilename,
		},
		ReportFile: FileReport{
			Path: config.ReportFilename,
		},
		ManifestFile: FileReport{
			Path:      config.ManifestFilename,
			SizeBytes: manifestFileSize,
			SizeHuman: bytefmt.ByteSize(uint64(manifestFileSize)),
		},
		PairedOutputs:       pairedReports,
		Profile:             ProfileReport{Directory: config.ProfileDir, CPUPath: config.CPUProfilePath, MemPath: config.MemProfilePath},
		Bucket:              bucketReport,
		Reads:               runStats.Reads,
		UncompressedBytes:   runStats.Bytes,
		OutputSizeBytes:     outputFileSize,
		SizeDifferenceBytes: sizeDifference,
		CompressionRatio:    compressionRatio,
		SizeReductionRatio:  sizeReductionRatio,
	}
	if err := WriteReport(report, config.ReportFilename); err != nil {
		return Result{}, err
	}

	slog.Debug("size reduced", "bytes", bytefmt.ByteSize(uint64(sizeDifference)), "ratio", sizeReductionRatio, "duration", timeDuration)
	return Result{Report: report}, nil
}
