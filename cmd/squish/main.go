package main

import (
	"context"
	"flag"
	"fmt"
	"log/slog"
	"os"
	"path/filepath"
	"squish"
)

func main() {
	squish.ConfigureLogging()

	_, sortMethodOptionStr := squish.GetSortingMethods()

	printVersion := flag.Bool("v", false, "print version information")
	sortMethodArg := flag.String("m", squish.DefaultSortMethod, "Fastq read sorting method. "+sortMethodOptionStr)
	cpuProfileFilename := flag.String("cpuProf", squish.DefaultCPUProfileFilename, "CPU profile filename")
	memProfileFilename := flag.String("memProf", squish.DefaultMemProfileFilename, "Memory profile filename")
	orderFilename := flag.String("orderFile", squish.DefaultOrderFilename, "File to record the order of sorted fastq reads")
	reportFilename := flag.String("reportFile", squish.DefaultReportFilename, "JSON report filename")
	manifestFilename := flag.String("manifestFile", squish.DefaultManifestFilename, "Text manifest filename listing absolute paths to output FASTQ files")
	outputDirArg := flag.String("outdir", squish.DefaultOutputDirNameBase, "Output dir")
	sortEngine := flag.String("engine", squish.DefaultSortEngine, "Sort engine. Options: memory, external")
	bucketStrategy := flag.String("bucket", squish.DefaultBucketStrategy, "External bucket strategy. Options: auto, sequence-prefix, quality-prefix, gc-range, hash, clump-minimizer")
	bucketCount := flag.Int("buckets", squish.DefaultExternalBucketCount, "External bucket count for bucket strategies that use a configurable count")
	clumpKmerLen := flag.Int("clumpK", squish.DefaultClumpKmerLen, "K-mer length used by the clump minimizer")
	tempDirArg := flag.String("tempdir", "tmp", "External bucket temp directory under the output dir")
	pairedFastqArg := flag.String("paired", "", "Comma- or semicolon-separated companion FASTQ files to reorder using the primary order file")
	pairedOutArg := flag.String("pairedOut", "", "Comma- or semicolon-separated output filenames for paired FASTQs under the output dir")
	checkPairs := flag.Bool("checkPairs", true, "Check companion FASTQ read names against the primary FASTQ before reordering")
	flag.Parse()

	if *printVersion {
		fmt.Println(squish.Version)
		return
	}

	cliArgs := flag.Args()
	if len(cliArgs) < 2 {
		slog.Error("not enough cli args provided", "required", 2, "got", len(cliArgs))
		flag.Usage()
		os.Exit(2)
	}

	config, err := configFromFlags(
		cliArgs,
		filepath.Clean(*outputDirArg),
		*sortMethodArg,
		*sortEngine,
		*bucketStrategy,
		*bucketCount,
		*clumpKmerLen,
		*orderFilename,
		*reportFilename,
		*manifestFilename,
		*cpuProfileFilename,
		*memProfileFilename,
		*tempDirArg,
		*pairedFastqArg,
		*pairedOutArg,
		*checkPairs,
	)
	if err != nil {
		slog.Error("invalid configuration", "error", err)
		os.Exit(2)
	}

	if _, err := squish.Run(context.Background(), config); err != nil {
		slog.Error("squish failed", "error", err)
		os.Exit(1)
	}
}

func configFromFlags(
	cliArgs []string,
	outputDir string,
	sortMethod string,
	sortEngine string,
	bucketStrategy string,
	bucketCount int,
	clumpKmerLen int,
	orderFilename string,
	reportFilename string,
	manifestFilename string,
	cpuProfileFilename string,
	memProfileFilename string,
	tempDirArg string,
	pairedFastqArg string,
	pairedOutArg string,
	checkPairs bool,
) (squish.Config, error) {
	inputFilepath := cliArgs[0]
	outputFilenameArg := cliArgs[1]

	outputFilepath, err := squish.OutputPath(outputDir, outputFilenameArg)
	if err != nil {
		return squish.Config{}, err
	}
	orderFilepath, err := squish.OutputPath(outputDir, orderFilename)
	if err != nil {
		return squish.Config{}, err
	}
	reportFilepath, err := squish.OutputPath(outputDir, reportFilename)
	if err != nil {
		return squish.Config{}, err
	}
	manifestFilepath, err := squish.OutputPath(outputDir, manifestFilename)
	if err != nil {
		return squish.Config{}, err
	}
	tempDir, err := squish.OutputPath(outputDir, tempDirArg)
	if err != nil {
		return squish.Config{}, err
	}

	pairedInputArgs := squish.SplitPathList(pairedFastqArg)
	pairedOutputArgs := squish.SplitPathList(pairedOutArg)
	if len(pairedOutputArgs) > 0 && len(pairedOutputArgs) != len(pairedInputArgs) {
		return squish.Config{}, fmt.Errorf("pairedOut has %d items but paired has %d items", len(pairedOutputArgs), len(pairedInputArgs))
	}
	if len(pairedOutputArgs) == 0 {
		for _, pairedInputArg := range pairedInputArgs {
			pairedOutputArgs = append(pairedOutputArgs, squish.DerivePairedOutputArg(pairedInputArg))
		}
	}
	pairedOutputPaths := make([]string, len(pairedOutputArgs))
	for i, pairedOutputArg := range pairedOutputArgs {
		pairedOutputPath, err := squish.OutputPath(outputDir, pairedOutputArg)
		if err != nil {
			return squish.Config{}, err
		}
		pairedOutputPaths[i] = pairedOutputPath
	}

	profileDir, err := squish.OutputPath(outputDir, squish.DefaultProfileDirnameBase+"."+sortMethod)
	if err != nil {
		return squish.Config{}, err
	}
	cpuProfilePath, err := squish.OutputPath(profileDir, cpuProfileFilename)
	if err != nil {
		return squish.Config{}, err
	}
	memProfilePath, err := squish.OutputPath(profileDir, memProfileFilename)
	if err != nil {
		return squish.Config{}, err
	}

	return squish.Config{
		SortMethod:            sortMethod,
		InputFilepath:         inputFilepath,
		OutputFilenameArg:     outputFilenameArg,
		OutputFilepath:        outputFilepath,
		PairedInputArgs:       pairedInputArgs,
		PairedInputFilepaths:  pairedInputArgs,
		PairedOutputArgs:      pairedOutputArgs,
		PairedOutputFilepaths: pairedOutputPaths,
		CheckPairs:            checkPairs,
		RecordDelim:           squish.RecordDelim,
		RecordHeaderChar:      squish.FastqHeaderChar,
		OrderFilename:         orderFilepath,
		ReportFilename:        reportFilepath,
		ManifestFilename:      manifestFilepath,
		OutputDir:             outputDir,
		SortEngine:            sortEngine,
		BucketStrategy:        bucketStrategy,
		BucketCount:           bucketCount,
		ClumpKmerLen:          clumpKmerLen,
		TempDir:               tempDir,
		ProfileDir:            profileDir,
		CPUProfilePath:        cpuProfilePath,
		MemProfilePath:        memProfilePath,
	}, nil
}
