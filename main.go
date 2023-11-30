package main

import (
	"code.cloudfoundry.org/bytefmt"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	_sort "squish/sort"
	"time"
)

// overwrite this at build time ;
// -ldflags="-X 'main.Version=someversion'"
var Version = "foo-version"

const fastqHeaderChar byte = '@'
const delim byte = '\n'

const defaultCpuProfileFilename string = "cpu.prof"
const defaultMemProfileFilename string = "mem.prof"
const defaultOrderFilename string = "order.txt"

// start CPU and Memory profiling
func Profiling(cpuFilename string, memFilename string) (*os.File, *os.File) {
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

// print the file size to the console log
func LogFileSize(filepath string, filetype string) int64 {
	inputFileSize, err := GetFileSize(filepath)
	if err != nil {
		log.Printf("WARNING: could not get size for file %v\n", filepath)
	}
	inputFileSizeBytes := bytefmt.ByteSize(uint64(inputFileSize))
	log.Printf("%v file %v of size %v Bytes\n", filetype, filepath, inputFileSizeBytes)
	return inputFileSize
}

// in case -v / --version CLI arg was used
func PrintVersionAndQuit() {
	fmt.Println(Version)
	os.Exit(0)
}

// make sure enough positional args were used
func MinCliPosArgs(args []string, n int) {
	if len(args) < n {
		log.Fatalf("Not enough cli args provided, %v args required, or use -h for help\n", n)
	}
}


func main() {
	timeStart := time.Now()
	sortMethodMap, sortMethodOptionStr := _sort.GetSortingMethods()
	defaultSortMethod, _ := _sort.GetDefaultSortMethod()

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
	config := _sort.Config{
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

	// run the sort methods
	config.Run()

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
