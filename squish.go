package squish

import (
	"context"
	"encoding/json"
	"fmt"
	"log/slog"
	"os"
	"path/filepath"
	"runtime/pprof"
	fastq "squish/fastq"
	_io "squish/fastqio"
	_sort "squish/sort"
	"strings"
	"time"

	"code.cloudfoundry.org/bytefmt"
)

// Version can be overwritten at build time with:
//
//	-ldflags="-X 'squish.Version=someversion'"
var Version = "foo-version"

const FastqHeaderChar byte = '@'
const RecordDelim byte = '\n'
const DefaultSortMethod = "alpha"
const DefaultSortDescription = "Alphabetical sort on sequence"
const DefaultCPUProfileFilename = "cpu.prof"
const DefaultMemProfileFilename = "mem.prof"
const DefaultOrderFilename = "order.txt"
const DefaultReportFilename = "report.json"
const DefaultManifestFilename = "manifest.txt"
const DefaultProfileDirnameBase = "profile"
const DefaultOutputDirNameBase = "output"
const DefaultSortEngine = "external"
const DefaultBucketStrategy = "auto"
const DefaultExternalBucketCount = 512
const DefaultClumpKmerLen = _sort.DefaultClumpKmerLen

type Result struct {
	Report Report
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

type SortDefinition struct {
	CLIArg      string
	Description string
	// Func is the legacy in-memory implementation.
	Func func(*[]fastq.FastqRead)
	// Strategy is the comparator shared by in-memory and external engines.
	Strategy _sort.SortStrategy
}

type Config struct {
	SortMethod            string
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
	ManifestFilename      string
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

func DefaultLogLevel() slog.Level {
	return slog.LevelDebug
}

func ConfigureLogging() {
	logger := slog.New(slog.NewTextHandler(os.Stderr, &slog.HandlerOptions{
		Level: DefaultLogLevel(),
	}))
	slog.SetDefault(logger)
}

func GetSortingMethods() (map[string]SortDefinition, string) {
	sortMethodMap := map[string]SortDefinition{
		DefaultSortMethod: SortDefinition{DefaultSortMethod, DefaultSortDescription, _sort.SortReadsSequence, _sort.AlphaSort{}},
		"gc":              SortDefinition{"gc", "GC Content Sort", _sort.SortReadsGC, _sort.GCSort{}},
		"qual":            SortDefinition{"qual", "Quality score sort", _sort.SortReadsQual, _sort.QualitySort{}},
		"alpha-heap":      SortDefinition{"alpha-heap", "Sequence alpha heap sort", _sort.HeapSortSequence, _sort.AlphaSort{}},
		"clump":           SortDefinition{"clump", "Clump-style read clustering for better gzip compression", _sort.SortReadsClump, _sort.ClumpSort{}},
	}

	sortMethodsDescr := map[string]string{}
	for key, value := range sortMethodMap {
		sortMethodsDescr[key] = value.Description
	}
	return sortMethodMap, fmt.Sprintf("Options: %v", sortMethodsDescr)
}

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

func normalizeConfig(config Config) (Config, SortDefinition, error) {
	if config.SortMethod == "" {
		config.SortMethod = DefaultSortMethod
	}
	if config.SortEngine == "" {
		config.SortEngine = DefaultSortEngine
	}
	if config.BucketStrategy == "" {
		config.BucketStrategy = DefaultBucketStrategy
	}
	if config.BucketCount == 0 {
		config.BucketCount = DefaultExternalBucketCount
	}
	if config.ClumpKmerLen == 0 {
		config.ClumpKmerLen = DefaultClumpKmerLen
	}
	if config.RecordDelim == 0 {
		config.RecordDelim = RecordDelim
	}
	if config.RecordHeaderChar == 0 {
		config.RecordHeaderChar = FastqHeaderChar
	}
	if config.TimeStart.IsZero() {
		config.TimeStart = time.Now()
	}
	if config.OutputDir == "" {
		config.OutputDir = DefaultOutputDirNameBase
	}
	config.OutputDir = filepath.Clean(config.OutputDir)

	sortMethodMap, _ := GetSortingMethods()
	sortDefinition, ok := sortMethodMap[config.SortMethod]
	if !ok {
		return Config{}, SortDefinition{}, fmt.Errorf("unknown sort method: %s", config.SortMethod)
	}
	if config.ClumpKmerLen < 1 {
		return Config{}, SortDefinition{}, fmt.Errorf("clumpK must be >= 1, got %d", config.ClumpKmerLen)
	}
	if sortDefinition.CLIArg == "clump" {
		sortDefinition.Func = func(reads *[]fastq.FastqRead) {
			_sort.SortReadsClumpK(reads, config.ClumpKmerLen)
		}
		sortDefinition.Strategy = _sort.ClumpSort{K: config.ClumpKmerLen}
	}

	if config.OutputFilenameArg == "" && config.OutputFilepath != "" {
		config.OutputFilenameArg = filepath.Base(config.OutputFilepath)
	}
	if config.OutputFilepath == "" {
		if config.OutputFilenameArg == "" {
			return Config{}, SortDefinition{}, fmt.Errorf("output filepath or output filename argument is required")
		}
		outputPath, err := OutputPath(config.OutputDir, config.OutputFilenameArg)
		if err != nil {
			return Config{}, SortDefinition{}, err
		}
		config.OutputFilepath = outputPath
	}
	if config.OrderFilename == "" {
		orderFilename, err := OutputPath(config.OutputDir, DefaultOrderFilename)
		if err != nil {
			return Config{}, SortDefinition{}, err
		}
		config.OrderFilename = orderFilename
	}
	if config.ReportFilename == "" {
		reportFilename, err := OutputPath(config.OutputDir, DefaultReportFilename)
		if err != nil {
			return Config{}, SortDefinition{}, err
		}
		config.ReportFilename = reportFilename
	}
	if config.ManifestFilename == "" {
		manifestFilename, err := OutputPath(config.OutputDir, DefaultManifestFilename)
		if err != nil {
			return Config{}, SortDefinition{}, err
		}
		config.ManifestFilename = manifestFilename
	}
	if config.TempDir == "" {
		tempDir, err := OutputPath(config.OutputDir, "tmp")
		if err != nil {
			return Config{}, SortDefinition{}, err
		}
		config.TempDir = tempDir
	}

	if len(config.PairedOutputArgs) > 0 && len(config.PairedOutputArgs) != len(config.PairedInputFilepaths) {
		return Config{}, SortDefinition{}, fmt.Errorf("paired output count %d does not match paired input count %d", len(config.PairedOutputArgs), len(config.PairedInputFilepaths))
	}
	if len(config.PairedOutputArgs) == 0 {
		for _, pairedInputArg := range config.PairedInputFilepaths {
			config.PairedOutputArgs = append(config.PairedOutputArgs, DerivePairedOutputArg(pairedInputArg))
		}
	}
	if len(config.PairedInputArgs) == 0 {
		config.PairedInputArgs = append(config.PairedInputArgs, config.PairedInputFilepaths...)
	}
	if len(config.PairedOutputFilepaths) == 0 {
		for _, pairedOutputArg := range config.PairedOutputArgs {
			pairedOutputPath, err := OutputPath(config.OutputDir, pairedOutputArg)
			if err != nil {
				return Config{}, SortDefinition{}, err
			}
			config.PairedOutputFilepaths = append(config.PairedOutputFilepaths, pairedOutputPath)
		}
	}
	if len(config.PairedOutputFilepaths) != len(config.PairedInputFilepaths) {
		return Config{}, SortDefinition{}, fmt.Errorf("paired output path count %d does not match paired input count %d", len(config.PairedOutputFilepaths), len(config.PairedInputFilepaths))
	}

	profileDir := config.ProfileDir
	if profileDir == "" {
		var err error
		profileDir, err = OutputPath(config.OutputDir, DefaultProfileDirnameBase+"."+sortDefinition.CLIArg)
		if err != nil {
			return Config{}, SortDefinition{}, err
		}
	}
	config.ProfileDir = profileDir
	if config.CPUProfilePath == "" {
		cpuProfilePath, err := OutputPath(profileDir, DefaultCPUProfileFilename)
		if err != nil {
			return Config{}, SortDefinition{}, err
		}
		config.CPUProfilePath = cpuProfilePath
	}
	if config.MemProfilePath == "" {
		memProfilePath, err := OutputPath(profileDir, DefaultMemProfileFilename)
		if err != nil {
			return Config{}, SortDefinition{}, err
		}
		config.MemProfilePath = memProfilePath
	}

	return config, sortDefinition, nil
}

func ensureOutputDirs(config Config) error {
	dirs := []string{
		config.OutputDir,
		filepath.Dir(config.OutputFilepath),
		filepath.Dir(config.OrderFilename),
		filepath.Dir(config.ReportFilename),
		filepath.Dir(config.ManifestFilename),
		config.ProfileDir,
	}
	for _, pairedOutputPath := range config.PairedOutputFilepaths {
		dirs = append(dirs, filepath.Dir(pairedOutputPath))
	}
	for _, dir := range dirs {
		if err := os.MkdirAll(dir, 0755); err != nil {
			return fmt.Errorf("create directory %q: %w", dir, err)
		}
	}
	return nil
}

func startProfiling(cpuFilename string, memFilename string) (*os.File, *os.File, error) {
	cpuFile, err := os.Create(cpuFilename)
	if err != nil {
		return nil, nil, fmt.Errorf("create CPU profile %q: %w", cpuFilename, err)
	}
	slog.Debug("CPU profile will be saved", "path", cpuFilename)
	if err := pprof.StartCPUProfile(cpuFile); err != nil {
		cpuFile.Close()
		return nil, nil, fmt.Errorf("start CPU profile: %w", err)
	}

	memFile, err := os.Create(memFilename)
	if err != nil {
		pprof.StopCPUProfile()
		cpuFile.Close()
		return nil, nil, fmt.Errorf("create memory profile %q: %w", memFilename, err)
	}
	slog.Debug("Memory profile will be saved", "path", memFilename)
	return cpuFile, memFile, nil
}

func stopProfiling(cpuFile *os.File, memFile *os.File) error {
	pprof.StopCPUProfile()
	var errs []string
	if err := pprof.WriteHeapProfile(memFile); err != nil {
		errs = append(errs, fmt.Sprintf("write heap profile: %v", err))
	}
	if err := cpuFile.Close(); err != nil {
		errs = append(errs, fmt.Sprintf("close CPU profile: %v", err))
	}
	if err := memFile.Close(); err != nil {
		errs = append(errs, fmt.Sprintf("close memory profile: %v", err))
	}
	if len(errs) > 0 {
		return fmt.Errorf(strings.Join(errs, "; "))
	}
	return nil
}

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
