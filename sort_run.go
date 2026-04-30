package squish

import (
	"fmt"
	"log/slog"
	"path/filepath"
	fastq "squish/fastq"
	_io "squish/fastqio"
	_sort "squish/sort"

	"code.cloudfoundry.org/bytefmt"
)

func RunSort(config Config, sortDefinition SortDefinition) (RunStats, error) {
	reader, err := _io.OpenReader(config.InputFilepath)
	if err != nil {
		return RunStats{}, err
	}
	defer reader.Close()

	writer, err := _io.OpenWriter(config.OutputFilepath)
	if err != nil {
		return RunStats{}, err
	}
	defer writer.Close()

	reads := []fastq.FastqRead{}
	totalByteSize, err := fastq.LoadReadsE(&reads, reader, &config.RecordDelim)
	if err != nil {
		return RunStats{}, err
	}
	slog.Info("reads loaded", "count", len(reads), "size", bytefmt.ByteSize(uint64(totalByteSize)))

	slog.Debug("starting read sort")
	sortDefinition.Func(&reads)
	slog.Debug("reads after sorting", "count", len(reads))

	slog.Debug("writing to output file", "path", config.OutputFilepath)
	if err := fastq.WriteReadsE(&reads, writer); err != nil {
		return RunStats{}, err
	}
	if err := fastq.SaveOrderE(&reads, config.OrderFilename); err != nil {
		return RunStats{}, err
	}

	return RunStats{Reads: len(reads), Bytes: totalByteSize}, nil
}

func RunPairedReorders(config Config, expectedReads int) ([]PairedRunStats, error) {
	if len(config.PairedInputFilepaths) == 0 {
		return nil, nil
	}

	var referenceNames []string
	var err error
	if config.CheckPairs {
		referenceNames, err = fastq.LoadNormalizedReadNames(config.InputFilepath, config.RecordDelim)
		if err != nil {
			return nil, fmt.Errorf("load primary read names for pairing checks: %w", err)
		}
		if len(referenceNames) != expectedReads {
			return nil, fmt.Errorf("primary read name count %d does not match sorted read count %d", len(referenceNames), expectedReads)
		}
	}

	pairedStats := make([]PairedRunStats, 0, len(config.PairedInputFilepaths))
	for i, inputPath := range config.PairedInputFilepaths {
		outputPath := config.PairedOutputFilepaths[i]
		slog.Debug("reordering paired fastq", "input", inputPath, "output", outputPath, "order", config.OrderFilename, "check_pairs", config.CheckPairs)

		stats, err := fastq.ReorderReadsByOrder(inputPath, outputPath, config.OrderFilename, config.RecordDelim, expectedReads, referenceNames)
		if err != nil {
			return nil, fmt.Errorf("reorder paired FASTQ %q: %w", inputPath, err)
		}
		inputSize := LogFileSize(inputPath, "Paired input")
		outputSize := LogFileSize(outputPath, "Paired output")
		pairedStats = append(pairedStats, PairedRunStats{
			InputArgument:     config.PairedInputArgs[i],
			InputPath:         inputPath,
			InputSizeBytes:    inputSize,
			OutputArgument:    config.PairedOutputArgs[i],
			OutputPath:        outputPath,
			Reads:             stats.Reads,
			UncompressedBytes: stats.Bytes,
			OutputSizeBytes:   outputSize,
		})
	}
	return pairedStats, nil
}

func RunExternalSort(config Config, sortDefinition SortDefinition) (RunStats, error) {
	tempDir := filepath.Join(config.TempDir, sortDefinition.CLIArg)
	sortConfig := _sort.ExternalBucketConfig{
		InputFilepath:  config.InputFilepath,
		OutputFilepath: config.OutputFilepath,
		OrderFilepath:  config.OrderFilename,
		TempDir:        tempDir,
		RecordDelim:    config.RecordDelim,
	}
	bucketer, err := GetBucketStrategy(config, sortDefinition)
	if err != nil {
		return RunStats{}, err
	}
	slog.Debug("starting external bucket sort", "sorter", sortDefinition.Strategy.Name(), "bucketer", bucketer.Name(), "buckets", bucketer.BucketCount(), "temp_dir", tempDir)
	stats, err := _sort.RunExternalBucketSort(sortConfig, sortDefinition.Strategy, bucketer)
	if err != nil {
		return RunStats{}, fmt.Errorf("external bucket sort failed: %w", err)
	}
	return RunStats{
		Reads:         stats.Reads,
		Bytes:         stats.Bytes,
		BucketsUsed:   stats.BucketsUsed,
		BucketCount:   stats.BucketCount,
		BucketName:    stats.BucketerName,
		BucketTempDir: stats.TempDir,
	}, nil
}

func GetBucketStrategy(config Config, sortDefinition SortDefinition) (_sort.BucketStrategy, error) {
	switch config.BucketStrategy {
	case "auto":
		return _sort.DefaultBucketStrategy(sortDefinition.Strategy, config.BucketCount), nil
	case "sequence-prefix":
		return _sort.NewSequencePrefixBuckets(2), nil
	case "quality-prefix":
		return _sort.NewQualityPrefixBuckets(1), nil
	case "gc-range":
		return _sort.NewGCRangeBuckets(config.BucketCount), nil
	case "hash":
		return _sort.NewHashBuckets(config.BucketCount), nil
	case "clump-minimizer":
		return _sort.NewClumpBuckets(config.BucketCount, config.ClumpKmerLen), nil
	default:
		return nil, fmt.Errorf("unknown bucket strategy: %s", config.BucketStrategy)
	}
}
