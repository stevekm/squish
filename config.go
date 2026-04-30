package squish

import (
	"fmt"
	"log/slog"
	"os"
	"path/filepath"
	fastq "squish/fastq"
	_sort "squish/sort"
	"time"
)

// Version can be overwritten at build time with:
//
//	-ldflags="-X 'squish.Version=someversion'"
var Version = "foo-version"

const FastqHeaderChar byte = '@'
const RecordDelim byte = '\n'
const DefaultSortMethod = "clump"
const DefaultSortDescription = "Clump-style read clustering for better gzip compression"
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
	ClumpMinCount         int  // filter pivot k-mers appearing fewer than this many times (0 = disabled)
	ClumpRComp            bool // reverse-complement minus-strand reads after clump sort
	ClumpRawPivot         bool // use lex-max canonical k-mer instead of max-hash pivot
	ClumpBorder           int  // bases excluded from each read end during pivot selection (Clumpify default: 1)
	QuantizeQuality       bool // bin quality scores to 4 Illumina levels after sorting (lossy)
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
		"alpha": SortDefinition{"alpha", "Alphabetical sort on sequence", _sort.SortReadsSequence, _sort.AlphaSort{}},
		"gc":    SortDefinition{"gc", "GC Content Sort", _sort.SortReadsGC, _sort.GCSort{}},
		"qual":  SortDefinition{"qual", "Quality score sort", _sort.SortReadsQual, _sort.QualitySort{}},
		"clump": SortDefinition{"clump", "Clump-style read clustering for better gzip compression", _sort.SortReadsClump, _sort.ClumpSort{}},
	}

	sortMethodsDescr := map[string]string{}
	for key, value := range sortMethodMap {
		sortMethodsDescr[key] = value.Description
	}
	return sortMethodMap, fmt.Sprintf("Options: %v", sortMethodsDescr)
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
		opts := _sort.ClumpSortOptions{
			K:        config.ClumpKmerLen,
			MinCount: config.ClumpMinCount,
			RComp:    config.ClumpRComp,
			RawPivot: config.ClumpRawPivot,
			Border:   config.ClumpBorder,
		}
		sortDefinition.Func = func(reads *[]fastq.FastqRead) {
			_sort.SortReadsClumpOpts(reads, opts)
		}
		sortDefinition.Strategy = _sort.ClumpSort{
			K:        config.ClumpKmerLen,
			MinCount: config.ClumpMinCount,
			RComp:    config.ClumpRComp,
			RawPivot: config.ClumpRawPivot,
			Border:   config.ClumpBorder,
		}
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
