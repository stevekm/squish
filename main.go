package main

import (
	"code.cloudfoundry.org/bytefmt"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	fastq "squish/fastq"
	_io "squish/io"
	_sort "squish/sort"
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
	log.Printf("CPU profile will be saved to %v\n", cpuFilename)
	pprof.StartCPUProfile(cpuFile)

	// Start memory profiling file
	memFile, err := os.Create(memFilename)
	if err != nil {
		log.Fatal(err)
	}
	log.Printf("Memory profile will be saved to %v\n", memFilename)

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
		log.Printf("WARNING: could not get size for file %v\n", filepath)
	}
	inputFileSizeBytes := bytefmt.ByteSize(uint64(inputFileSize))
	log.Printf("%v file %v of size %v Bytes\n", filetype, filepath, inputFileSizeBytes)
	return inputFileSize
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

func RunSort(config Config) {
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
	fastq.LoadReads(&reads, reader, &config.RecordDelim)
	log.Printf("%v reads loaded\n", len(reads))

	// sort the fastq reads
	log.Printf("starting read sort")
	config.SortMethod.Func(&reads)
	log.Printf("%v reads after sorting\n", len(reads))

	// write the fastq reads
	log.Printf("Writing to output file %v\n", config.OutputFilepath)
	fastq.WriteReads(&reads, writer)

	// save the order of the sorted reads to file
	fastq.SaveOrder(&reads, config.OrderFilename)
}

func GetSortingMethods() (map[string]SortMethod, string) {
	// parse the available sorting methods
	sortMethodMap := map[string]SortMethod{
		defaultSortMethod: SortMethod{defaultSortMethod, defaultSortDescrption, _sort.SortReadsSequence}, // alpha
		"gc":              SortMethod{"gc", "GC Content Sort", _sort.SortReadsGC},
		"qual":            SortMethod{"qual", "Quality score sort", _sort.SortReadsQual},
		"alpha-heap":      SortMethod{"alpha-heap", "Sequence alpha heap sort", _sort.HeapSortSequence},
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
	Func        func(*[]fastq.FastqRead)
}

type Config struct {
	SortMethod       SortMethod
	InputFilepath    string
	InputFileSize    int64
	OutputFilepath   string
	RecordDelim      byte
	RecordHeaderChar byte
	TimeStart        time.Time
	OrderFilename    string
}

func main() {
	timeStart := time.Now()
	sortMethodMap, sortMethodOptionStr := GetSortingMethods()

	// get command line args
	printVersion := flag.Bool("v", false, "print version information")
	sortMethodArg := flag.String("m", defaultSortMethod, "Fastq read sorting method. "+sortMethodOptionStr)
	cpuProfileFilename := flag.String("cpuProf", defaultCpuProfileFilename, "CPU profile filename")
	memProfileFilename := flag.String("memProf", defaultMemProfileFilename, "Memory profile filename")
	orderFilename := flag.String("orderFile", defaultOrderFilename, "File to record the order of sorted fastq reads")
	flag.Parse()
	cliArgs := flag.Args() // all positional args passed
	if *printVersion {
		PrintVersionAndQuit()
	}
	MinCliPosArgs(cliArgs, 2)

	// get positional cli args
	inputFilepath := cliArgs[0]
	outputFilepath := cliArgs[1]

	inputFileSize := LogFileSize(inputFilepath, "Input")

	// initialize config
	_, ok := sortMethodMap[*sortMethodArg]
	if !ok {
		log.Fatalf("ERROR: Unknown sort method: %v\n", *sortMethodArg)
	}
	config := Config{
		SortMethod:       sortMethodMap[*sortMethodArg],
		InputFilepath:    inputFilepath,
		InputFileSize:    inputFileSize,
		OutputFilepath:   outputFilepath,
		RecordDelim:      delim,
		RecordHeaderChar: fastqHeaderChar,
		TimeStart:        timeStart,
		OrderFilename:    *orderFilename,
	}

	//
	//
	//
	// start profiler
	cpuFile, memFile := Profiling(*cpuProfileFilename, *memProfileFilename)
	defer cpuFile.Close()
	defer memFile.Close()
	defer pprof.StopCPUProfile()
	defer pprof.WriteHeapProfile(memFile)

	// insert sort methods here
	log.Printf("Using sort method: %v\n", config.SortMethod.CLIArg)
	RunSort(config)

	//
	//
	//

	// print some stuff to the console log
	timeStop := time.Now()
	timeDuration := timeStop.Sub(config.TimeStart)
	outputFileSize := LogFileSize(config.OutputFilepath, "Output")
	sizeDifference := config.InputFileSize - outputFileSize
	sizeDifferenceBytes := bytefmt.ByteSize(uint64(sizeDifference))

	log.Printf("Size reduced by %v Bytes (%.4f) in %v\n",
		sizeDifferenceBytes,
		float64(sizeDifference)/float64(inputFileSize),
		timeDuration,
	)
}
