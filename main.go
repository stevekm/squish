package main

import (
	"fmt"
	_io "squish/io"
	_sort "squish/sort"
	fastq "squish/fastq"
	"bufio"
	"code.cloudfoundry.org/bytefmt"
	"flag"
	"log"
	"os"
	"runtime/pprof"
	"strconv"
	"time"
)

// overwrite this at build time ;
// -ldflags="-X 'main.Version=someversion'"
var Version = "foo-version"

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

func SaveOrder(readsBuffer *[]fastq.FastqRead) {
	outputFilepath := "order.txt"
	log.Printf("Saving read order to file %v for %v reads\n", outputFilepath, len(*readsBuffer))
	outputFile, err := os.Create(outputFilepath)
	if err != nil {
		log.Fatalf("Error creating output file: %v\n", err)
	}
	defer outputFile.Close()

	writer := bufio.NewWriter(outputFile)
	// defer writer.Close()
	defer writer.Flush()
	for _, read := range *readsBuffer {
		_, err := writer.WriteString(strconv.Itoa(read.I) + "\n")
		if err != nil {
			log.Fatalf("Error writing to file: %v\n", err)
		}
	}
}

func PrintVersionAndQuit() {
	fmt.Println(Version)
	os.Exit(0)
}

func MinCliPosArgs(args []string, n int) {
	if len(args) < n {
		log.Fatalf("Not enough cli args provided, %v args required\n", n)
	}
}

func RunAlphaSort(config Config) {
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
	_sort.SortReads(&reads)
	log.Printf("%v reads after sorting\n", len(reads))

	// write the fastq reads
	log.Printf("Writing to output file %v\n", config.OutputFilepath)
	fastq.WriteReads(&reads, writer)

	// save the order of the sorted reads to file
	SaveOrder(&reads)
}

func RunGCSort(config Config) {
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
	_sort.SortReadsGC(&reads)
	log.Printf("%v reads after sorting\n", len(reads))

	// write the fastq reads
	log.Printf("Writing to output file %v\n", config.OutputFilepath)
	fastq.WriteReads(&reads, writer)

	// save the order of the sorted reads to file
	SaveOrder(&reads)
}

func RunQualSort(config Config) {
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
	_sort.SortReadsQual(&reads)
	log.Printf("%v reads after sorting\n", len(reads))

	// write the fastq reads
	log.Printf("Writing to output file %v\n", config.OutputFilepath)
	fastq.WriteReads(&reads, writer)

	// save the order of the sorted reads to file
	SaveOrder(&reads)
}

type Config struct {
	SortMethod       string
	SortMethods      []string
	InputFilepath    string
	InputFileSize    int64
	OutputFilepath   string
	RecordDelim      byte
	RecordHeaderChar byte
	TimeStart        time.Time
}

const fastqHeaderChar byte = '@'
const delim byte = '\n'
const defaultSortMethod string = "alpha"
const cpuProfileFilename string = "cpu.prof"
const memProfileFilename string = "mem.prof"

func main() {
	timeStart := time.Now()
	sortMethods := []string{defaultSortMethod, "gc", "qual"}
	sortMethodOptionStr := fmt.Sprintf("Options: %v", sortMethods)

	// start profiler
	cpuFile, memFile := Profiling(cpuProfileFilename, memProfileFilename)
	defer cpuFile.Close()
	defer memFile.Close()
	defer pprof.StopCPUProfile()
	defer pprof.WriteHeapProfile(memFile)

	// get command line args
	printVersion := flag.Bool("v", false, "print version information")
	sortMethod := flag.String("m", defaultSortMethod, "Fastq read sorting method. "+sortMethodOptionStr)
	flag.Parse()
	cliArgs := flag.Args() // all positional args passed
	if *printVersion {
		PrintVersionAndQuit()
	}
	MinCliPosArgs(cliArgs, 2)
	inputFilepath := cliArgs[0]
	outputFilepath := cliArgs[1]

	inputFileSize := LogFileSize(inputFilepath, "Input")
	config := Config{
		SortMethod:       *sortMethod,
		SortMethods:      sortMethods,
		InputFilepath:    inputFilepath,
		InputFileSize:    inputFileSize,
		OutputFilepath:   outputFilepath,
		RecordDelim:      delim,
		RecordHeaderChar: fastqHeaderChar,
		TimeStart:        timeStart,
	}

	//
	//
	//

	// insert modular sort methods here
	log.Printf("Using sort method: %v\n", config.SortMethod)
	switch config.SortMethod {
	case "alpha":
		RunAlphaSort(config)
	case "gc":
		RunGCSort(config)
	case "qual":
		RunQualSort(config)
	default:
		log.Printf("Using default sort method: %v\n", defaultSortMethod)
		RunAlphaSort(config)
	}

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
