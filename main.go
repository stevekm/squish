package main

import (
	"code.cloudfoundry.org/bytefmt"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"runtime/pprof"
	fastq "squish/fastq"
	_io "squish/io"
	_sort "squish/sort"
	"strings"
	"time"
)

// overwrite this at build time ;
// -ldflags="-X 'main.Version=someversion'"
var Version = "foo-version"

const fastqHeaderChar byte = '@'
const delim byte = '\n'
const defaultSortMethod string = "alpha"
const defaultSortDescrption string = "Alphabetical sort on sequence"

// const defaultSortFunc func() = _sort.SortReadsSequence
const defaultCpuProfileFilename string = "cpu.prof"
const defaultMemProfileFilename string = "mem.prof"
const defaultOrderFilename string = "order.txt"
const defaultReportFilename string = "report.json"
const defaultProfileDirnameBase string = "profile"
const defaultOutputDirNameBase string = "output"
const defaultSortEngine string = "memory"
const defaultBucketStrategy string = "auto"
const defaultExternalBucketCount int = 512
const defaultClumpKmerLen int = _sort.DefaultClumpKmerLen

func DefaultLogLevel() slog.Level {
	return slog.LevelDebug
}

func ConfigureLogging() {
	logger := slog.New(slog.NewTextHandler(os.Stderr, &slog.HandlerOptions{
		Level: DefaultLogLevel(),
	}))
	slog.SetDefault(logger)
}

func Profiling(cpuFilename string, memFilename string) (*os.File, *os.File) {
	// start CPU and Memory profiling
	//
	// NOTES:
	// https://pkg.go.dev/runtime/pprof
	// https://github.com/google/pprof/blob/main/doc/README.md
	// $ go tool pprof cpu.prof
	// $ go tool pprof mem.prof
	// (pprof) top
	//
	// caller must also run these;
	// defer cpuFile.Close()
	// defer memFile.Close()
	// defer pprof.StopCPUProfile()
	// defer pprof.WriteHeapProfile(memFile)
	//
	// see also; https://pkg.go.dev/net/http/pprof

	// Start CPU profiling
	cpuFile, err := os.Create(cpuFilename)
	if err != nil {
		log.Fatal(err)
	}
	slog.Debug("CPU profile will be saved", "path", cpuFilename)
	pprof.StartCPUProfile(cpuFile)

	// Start memory profiling file
	memFile, err := os.Create(memFilename)
	if err != nil {
		log.Fatal(err)
	}
	slog.Debug("Memory profile will be saved", "path", memFilename)

	return cpuFile, memFile
}

func GetFileSize(filepath string) (int64, error) {
	fi, err := os.Stat(filepath)
	if err != nil {
		return 0, err
	}
	return fi.Size(), nil
}

func LogFileSize(filepath string, filetype string) int64 {
	inputFileSize, err := GetFileSize(filepath)
	if err != nil {
		slog.Debug("could not get size for file", "path", filepath, "error", err)
	}
	inputFileSizeBytes := bytefmt.ByteSize(uint64(inputFileSize))
	slog.Debug("file size", "type", filetype, "path", filepath, "size", inputFileSizeBytes)
	return inputFileSize
}

type RunStats struct {
	Reads         int
	Bytes         int
	BucketsUsed   int
	BucketCount   int
	BucketName    string
	BucketTempDir string
}

type PairedRunStats struct {
	InputArgument     string
	InputPath         string
	InputSizeBytes    int64
	OutputArgument    string
	OutputPath        string
	Reads             int
	UncompressedBytes int
	OutputSizeBytes   int64
}

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

func WriteReport(report Report, reportPath string) {
	reportJSON, err := json.MarshalIndent(report, "", "  ")
	if err != nil {
		log.Fatalf("ERROR: could not marshal report: %v\n", err)
	}
	if err := os.WriteFile(reportPath, reportJSON, 0644); err != nil {
		log.Fatalf("ERROR: could not write report file %v: %v\n", reportPath, err)
	}
	slog.Debug("report written", "path", reportPath)
}

func PrintVersionAndQuit() {
	fmt.Println(Version)
	os.Exit(0)
}

func MinCliPosArgs(args []string, n int) {
	if len(args) < n {
		log.Fatalf("Not enough cli args provided, %v args required, or use -h for help\n", n)
	}
}

func OutputPath(outputDir string, itemPath string) string {
	// All user-supplied output filenames are interpreted relative to outputDir.
	// The Rel check prevents paths such as "../outside.fastq.gz" from escaping
	// the run directory.
	outputPath := filepath.Join(outputDir, itemPath)
	relPath, err := filepath.Rel(outputDir, outputPath)
	if err != nil || relPath == ".." || strings.HasPrefix(relPath, ".."+string(os.PathSeparator)) {
		log.Fatalf("ERROR: output item path %v must stay within output directory %v\n", itemPath, outputDir)
	}
	return outputPath
}

func AbsolutePath(path string) string {
	absolutePath, err := filepath.Abs(path)
	if err != nil {
		log.Fatalf("ERROR: could not resolve absolute path for %v: %v\n", path, err)
	}
	return absolutePath
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

func RunSort(config Config) RunStats {
	// run the chosen sorting method on the fastq file

	// input
	reader := _io.GetReader(config.InputFilepath)
	defer reader.Close()

	// output
	writer := _io.GetWriter(config.OutputFilepath)
	defer writer.Close()

	// hold all the reads from the file in here
	reads := []fastq.FastqRead{}

	// load all reads from file
	totalByteSize := fastq.LoadReads(&reads, reader, &config.RecordDelim)
	slog.Info("reads loaded", "count", len(reads), "size", bytefmt.ByteSize(uint64(totalByteSize)))

	// sort the fastq reads
	slog.Debug("starting read sort")
	config.SortMethod.Func(&reads)
	slog.Debug("reads after sorting", "count", len(reads))

	// write the fastq reads
	slog.Debug("writing to output file", "path", config.OutputFilepath)
	fastq.WriteReads(&reads, writer)

	// save the order of the sorted reads to file
	fastq.SaveOrder(&reads, config.OrderFilename)

	return RunStats{
		Reads: len(reads),
		Bytes: totalByteSize,
	}
}

func RunPairedReorders(config Config, expectedReads int) []PairedRunStats {
	if len(config.PairedInputFilepaths) == 0 {
		return nil
	}

	var referenceNames []string
	var err error
	if config.CheckPairs {
		// Load normalized R1 names once and compare every companion FASTQ in
		// original input order before applying the sorted order file.
		referenceNames, err = fastq.LoadNormalizedReadNames(config.InputFilepath, config.RecordDelim)
		if err != nil {
			log.Fatalf("ERROR: could not load primary read names for pairing checks: %v\n", err)
		}
		if len(referenceNames) != expectedReads {
			log.Fatalf("ERROR: primary read name count %d does not match sorted read count %d\n", len(referenceNames), expectedReads)
		}
	}

	pairedStats := make([]PairedRunStats, 0, len(config.PairedInputFilepaths))
	for i, inputPath := range config.PairedInputFilepaths {
		outputPath := config.PairedOutputFilepaths[i]
		slog.Debug(
			"reordering paired fastq",
			"input", inputPath,
			"output", outputPath,
			"order", config.OrderFilename,
			"check_pairs", config.CheckPairs,
		)

		stats, err := fastq.ReorderReadsByOrder(inputPath, outputPath, config.OrderFilename, config.RecordDelim, expectedReads, referenceNames)
		if err != nil {
			log.Fatalf("ERROR: could not reorder paired FASTQ %v: %v\n", inputPath, err)
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
	return pairedStats
}

func RunExternalSort(config Config) RunStats {
	// Keep temporary buckets under the configured output directory and isolate
	// each sort method in its own temp subdirectory.
	tempDir := filepath.Join(config.TempDir, config.SortMethod.CLIArg)
	sortConfig := _sort.ExternalBucketConfig{
		InputFilepath:  config.InputFilepath,
		OutputFilepath: config.OutputFilepath,
		OrderFilepath:  config.OrderFilename,
		TempDir:        tempDir,
		RecordDelim:    config.RecordDelim,
	}
	// Bucket selection is deliberately separate from the sort strategy so the
	// CLI can compose different file-backed layouts with the same comparator.
	bucketer := GetBucketStrategy(config)
	slog.Debug(
		"starting external bucket sort",
		"sorter", config.SortMethod.Strategy.Name(),
		"bucketer", bucketer.Name(),
		"buckets", bucketer.BucketCount(),
		"temp_dir", tempDir,
	)
	stats, err := _sort.RunExternalBucketSort(sortConfig, config.SortMethod.Strategy, bucketer)
	if err != nil {
		log.Fatalf("ERROR: external bucket sort failed: %v\n", err)
	}
	return RunStats{
		Reads:         stats.Reads,
		Bytes:         stats.Bytes,
		BucketsUsed:   stats.BucketsUsed,
		BucketCount:   stats.BucketCount,
		BucketName:    stats.BucketerName,
		BucketTempDir: stats.TempDir,
	}
}

func GetSortingMethods() (map[string]SortMethod, string) {
	// parse the available sorting methods
	sortMethodMap := map[string]SortMethod{
		defaultSortMethod: SortMethod{defaultSortMethod, defaultSortDescrption, _sort.SortReadsSequence, _sort.AlphaSort{}}, // alpha
		"gc":              SortMethod{"gc", "GC Content Sort", _sort.SortReadsGC, _sort.GCSort{}},
		"qual":            SortMethod{"qual", "Quality score sort", _sort.SortReadsQual, _sort.QualitySort{}},
		"alpha-heap":      SortMethod{"alpha-heap", "Sequence alpha heap sort", _sort.HeapSortSequence, _sort.AlphaSort{}},
		"clump":           SortMethod{"clump", "Clump-style read clustering for better gzip compression", _sort.SortReadsClump, _sort.ClumpSort{}},
	}
	// minimal map for help text printing
	sortMethodsDescr := map[string]string{}
	for key, value := range sortMethodMap {
		sortMethodsDescr[key] = value.Description
	}
	// help text
	sortMethodOptionStr := fmt.Sprintf("Options: %v", sortMethodsDescr)
	return sortMethodMap, sortMethodOptionStr
}

type SortMethod struct {
	CLIArg      string
	Description string
	// Func is the legacy in-memory implementation.
	Func func(*[]fastq.FastqRead)
	// Strategy is the comparator shared by in-memory and external engines.
	Strategy _sort.SortStrategy
}

type Config struct {
	SortMethod            SortMethod
	InputFilepath         string
	InputFileSize         int64
	OutputFilenameArg     string
	OutputFilepath        string
	PairedInputArgs       []string
	PairedInputFilepaths  []string
	PairedOutputArgs      []string
	PairedOutputFilepaths []string
	CheckPairs            bool
	RecordDelim           byte
	RecordHeaderChar      byte
	TimeStart             time.Time
	OrderFilename         string
	ReportFilename        string
	OutputDir             string
	SortEngine            string
	BucketStrategy        string
	BucketCount           int
	ClumpKmerLen          int
	TempDir               string
	ProfileDir            string
	CPUProfilePath        string
	MemProfilePath        string
}

func GetBucketStrategy(config Config) _sort.BucketStrategy {
	// "auto" chooses the safest default for the selected sorter. Explicit
	// strategies are useful for experiments and profiling.
	switch config.BucketStrategy {
	case "auto":
		return _sort.DefaultBucketStrategy(config.SortMethod.Strategy, config.BucketCount)
	case "sequence-prefix":
		return _sort.NewSequencePrefixBuckets(2)
	case "quality-prefix":
		return _sort.NewQualityPrefixBuckets(1)
	case "gc-range":
		return _sort.NewGCRangeBuckets(config.BucketCount)
	case "hash":
		return _sort.NewHashBuckets(config.BucketCount)
	case "clump-minimizer":
		return _sort.NewClumpBuckets(config.BucketCount, config.ClumpKmerLen)
	default:
		log.Fatalf("ERROR: Unknown bucket strategy: %v\n", config.BucketStrategy)
	}
	return nil
}

func main() {
	ConfigureLogging()

	timeStart := time.Now()
	sortMethodMap, sortMethodOptionStr := GetSortingMethods()

	// get command line args
	printVersion := flag.Bool("v", false, "print version information")
	sortMethodArg := flag.String("m", defaultSortMethod, "Fastq read sorting method. "+sortMethodOptionStr)
	cpuProfileFilename := flag.String("cpuProf", defaultCpuProfileFilename, "CPU profile filename")
	memProfileFilename := flag.String("memProf", defaultMemProfileFilename, "Memory profile filename")
	orderFilename := flag.String("orderFile", defaultOrderFilename, "File to record the order of sorted fastq reads")
	reportFilename := flag.String("reportFile", defaultReportFilename, "JSON report filename")
	outputDirArg := flag.String("outdir", defaultOutputDirNameBase, "Output dir")
	sortEngine := flag.String("engine", defaultSortEngine, "Sort engine. Options: memory, external")
	bucketStrategy := flag.String("bucket", defaultBucketStrategy, "External bucket strategy. Options: auto, sequence-prefix, quality-prefix, gc-range, hash, clump-minimizer")
	bucketCount := flag.Int("buckets", defaultExternalBucketCount, "External bucket count for bucket strategies that use a configurable count")
	clumpKmerLen := flag.Int("clumpK", defaultClumpKmerLen, "K-mer length used by the clump minimizer")
	tempDirArg := flag.String("tempdir", "tmp", "External bucket temp directory under the output dir")
	pairedFastqArg := flag.String("paired", "", "Comma- or semicolon-separated companion FASTQ files to reorder using the primary order file")
	pairedOutArg := flag.String("pairedOut", "", "Comma- or semicolon-separated output filenames for paired FASTQs under the output dir")
	checkPairs := flag.Bool("checkPairs", true, "Check companion FASTQ read names against the primary FASTQ before reordering")
	flag.Parse()
	cliArgs := flag.Args() // all positional args passed
	if *printVersion {
		PrintVersionAndQuit()
	}
	MinCliPosArgs(cliArgs, 2)

	// get positional cli args
	inputFilepath := cliArgs[0]
	outputFilenameArg := cliArgs[1]
	outputDir := filepath.Clean(*outputDirArg)
	outputFilepath := OutputPath(outputDir, outputFilenameArg)
	orderFilepath := OutputPath(outputDir, *orderFilename)
	reportFilepath := OutputPath(outputDir, *reportFilename)
	tempDir := OutputPath(outputDir, *tempDirArg)
	pairedInputArgs := SplitPathList(*pairedFastqArg)
	pairedOutputArgs := SplitPathList(*pairedOutArg)
	if len(pairedOutputArgs) > 0 && len(pairedOutputArgs) != len(pairedInputArgs) {
		log.Fatalf("ERROR: pairedOut has %d items but paired has %d items\n", len(pairedOutputArgs), len(pairedInputArgs))
	}
	if len(pairedOutputArgs) == 0 {
		for _, pairedInputArg := range pairedInputArgs {
			pairedOutputArgs = append(pairedOutputArgs, DerivePairedOutputArg(pairedInputArg))
		}
	}
	pairedOutputPaths := make([]string, len(pairedOutputArgs))
	for i, pairedOutputArg := range pairedOutputArgs {
		pairedOutputPaths[i] = OutputPath(outputDir, pairedOutputArg)
	}

	inputFileSize := LogFileSize(inputFilepath, "Input")

	// initialize config
	sortMethod, ok := sortMethodMap[*sortMethodArg]
	if !ok {
		log.Fatalf("ERROR: Unknown sort method: %v\n", *sortMethodArg)
	}
	if *clumpKmerLen < 1 {
		log.Fatalf("ERROR: clumpK must be >= 1, got %v\n", *clumpKmerLen)
	}
	if sortMethod.CLIArg == "clump" {
		// The clump k-mer length is configured after flag parsing, so replace
		// the default clump function and strategy with k-aware versions here.
		sortMethod.Func = func(reads *[]fastq.FastqRead) {
			_sort.SortReadsClumpK(reads, *clumpKmerLen)
		}
		sortMethod.Strategy = _sort.ClumpSort{K: *clumpKmerLen}
	}
	config := Config{
		SortMethod:            sortMethod,
		InputFilepath:         inputFilepath,
		InputFileSize:         inputFileSize,
		OutputFilenameArg:     outputFilenameArg,
		OutputFilepath:        outputFilepath,
		PairedInputArgs:       pairedInputArgs,
		PairedInputFilepaths:  pairedInputArgs,
		PairedOutputArgs:      pairedOutputArgs,
		PairedOutputFilepaths: pairedOutputPaths,
		CheckPairs:            *checkPairs,
		RecordDelim:           delim,
		RecordHeaderChar:      fastqHeaderChar,
		TimeStart:             timeStart,
		OrderFilename:         orderFilepath,
		ReportFilename:        reportFilepath,
		OutputDir:             outputDir,
		SortEngine:            *sortEngine,
		BucketStrategy:        *bucketStrategy,
		BucketCount:           *bucketCount,
		ClumpKmerLen:          *clumpKmerLen,
		TempDir:               tempDir,
	}

	if err := os.MkdirAll(outputDir, 0755); err != nil {
		log.Fatalf("ERROR: could not create output directory %v: %v\n", outputDir, err)
	}
	slog.Debug("output directory", "path", outputDir)

	outputFileDir := filepath.Dir(config.OutputFilepath)
	if err := os.MkdirAll(outputFileDir, 0755); err != nil {
		log.Fatalf("ERROR: could not create output file directory %v: %v\n", outputFileDir, err)
	}

	orderFileDir := filepath.Dir(config.OrderFilename)
	if err := os.MkdirAll(orderFileDir, 0755); err != nil {
		log.Fatalf("ERROR: could not create order file directory %v: %v\n", orderFileDir, err)
	}

	reportFileDir := filepath.Dir(config.ReportFilename)
	if err := os.MkdirAll(reportFileDir, 0755); err != nil {
		log.Fatalf("ERROR: could not create report file directory %v: %v\n", reportFileDir, err)
	}
	for _, pairedOutputPath := range config.PairedOutputFilepaths {
		pairedOutputDir := filepath.Dir(pairedOutputPath)
		if err := os.MkdirAll(pairedOutputDir, 0755); err != nil {
			log.Fatalf("ERROR: could not create paired output file directory %v: %v\n", pairedOutputDir, err)
		}
	}

	profileDirNameBase := OutputPath(outputDir, defaultProfileDirnameBase+"."+config.SortMethod.CLIArg)
	if err := os.MkdirAll(profileDirNameBase, 0755); err != nil {
		log.Fatalf("ERROR: could not create profile directory %v: %v\n", profileDirNameBase, err)
	}
	slog.Debug("saving profile", "path", profileDirNameBase)

	cpuProfilePath := OutputPath(profileDirNameBase, *cpuProfileFilename)
	memProfilePath := OutputPath(profileDirNameBase, *memProfileFilename)
	config.ProfileDir = profileDirNameBase
	config.CPUProfilePath = cpuProfilePath
	config.MemProfilePath = memProfilePath

	//
	//
	//
	// start profiler
	cpuFile, memFile := Profiling(cpuProfilePath, memProfilePath)
	defer cpuFile.Close()
	defer memFile.Close()
	defer pprof.StopCPUProfile()
	defer pprof.WriteHeapProfile(memFile)

	// insert sort methods here
	slog.Debug("using sort method", "method", config.SortMethod.CLIArg, "clump_k", config.ClumpKmerLen)
	var runStats RunStats
	switch config.SortEngine {
	case "memory":
		runStats = RunSort(config)
	case "external":
		runStats = RunExternalSort(config)
	default:
		log.Fatalf("ERROR: Unknown sort engine: %v\n", config.SortEngine)
	}
	pairedStats := RunPairedReorders(config, runStats.Reads)

	//
	//
	//

	// print some stuff to the console log
	timeStop := time.Now()
	timeDuration := timeStop.Sub(config.TimeStart)
	outputFileSize := LogFileSize(config.OutputFilepath, "Output")
	sizeDifference := config.InputFileSize - outputFileSize
	sizeDifferenceBytes := bytefmt.ByteSize(uint64(sizeDifference))
	compressionRatio := 0.0
	sizeReductionRatio := 0.0
	if config.InputFileSize > 0 {
		compressionRatio = float64(outputFileSize) / float64(config.InputFileSize)
		sizeReductionRatio = float64(sizeDifference) / float64(config.InputFileSize)
	}
	var bucketReport *BucketReport
	if config.SortEngine == "external" {
		// Recreate the selected bucket strategy only for reporting whether the
		// emitted bucket order is a true global ordering for this sorter.
		bucketer := GetBucketStrategy(config)
		bucketReport = &BucketReport{
			Strategy:   runStats.BucketName,
			Count:      runStats.BucketCount,
			Used:       runStats.BucketsUsed,
			TempDir:    runStats.BucketTempDir,
			OrderedFor: bucketer.OrderedFor(config.SortMethod.Strategy),
		}
	}
	pairedReports := make([]PairedReport, 0, len(pairedStats))
	for _, pairedStat := range pairedStats {
		pairedReports = append(pairedReports, PairedReport{
			Input: FileReport{
				Path:      pairedStat.InputPath,
				Argument:  pairedStat.InputArgument,
				SizeBytes: pairedStat.InputSizeBytes,
				SizeHuman: bytefmt.ByteSize(uint64(pairedStat.InputSizeBytes)),
			},
			Output: FileReport{
				Path:      AbsolutePath(pairedStat.OutputPath),
				Argument:  pairedStat.OutputArgument,
				SizeBytes: pairedStat.OutputSizeBytes,
				SizeHuman: bytefmt.ByteSize(uint64(pairedStat.OutputSizeBytes)),
			},
			Reads:             pairedStat.Reads,
			UncompressedBytes: pairedStat.UncompressedBytes,
			OutputSizeBytes:   pairedStat.OutputSizeBytes,
		})
	}

	slog.Debug(
		"size reduced",
		"bytes", sizeDifferenceBytes,
		"ratio", sizeReductionRatio,
		"duration", timeDuration,
	)

	WriteReport(Report{
		Version:              Version,
		StartedAt:            config.TimeStart.Format(time.RFC3339Nano),
		FinishedAt:           timeStop.Format(time.RFC3339Nano),
		Duration:             timeDuration.String(),
		DurationMilliseconds: timeDuration.Milliseconds(),
		SortMethod:           config.SortMethod.CLIArg,
		SortDescription:      config.SortMethod.Description,
		SortEngine:           config.SortEngine,
		ClumpKmerLength:      config.ClumpKmerLen,
		Input: FileReport{
			Path:      config.InputFilepath,
			SizeBytes: config.InputFileSize,
			SizeHuman: bytefmt.ByteSize(uint64(config.InputFileSize)),
		},
		Output: FileReport{
			Path:      AbsolutePath(config.OutputFilepath),
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
		PairedOutputs: pairedReports,
		Profile: ProfileReport{
			Directory: config.ProfileDir,
			CPUPath:   config.CPUProfilePath,
			MemPath:   config.MemProfilePath,
		},
		Bucket:              bucketReport,
		Reads:               runStats.Reads,
		UncompressedBytes:   runStats.Bytes,
		OutputSizeBytes:     outputFileSize,
		SizeDifferenceBytes: sizeDifference,
		CompressionRatio:    compressionRatio,
		SizeReductionRatio:  sizeReductionRatio,
	}, config.ReportFilename)
}
