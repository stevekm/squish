package main

import (
	"fmt"
	// "compress/gzip"
	"bufio"
	"code.cloudfoundry.org/bytefmt"
	"flag"
	gzip "github.com/klauspost/pgzip"
	"log"
	"os"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"time"
)

// overwrite this at build time ;
// -ldflags="-X 'main.Version=someversion'"
var Version = "foo-version"

// NOTE: fastq read methods https://pkg.go.dev/github.com/biogo/biogo/io/seqio/fastq#Reader.Read
// https://en.wikipedia.org/wiki/FASTQ_format
type FastqRead struct {
	Id            []byte
	Sequence      []byte
	Plus          []byte
	QualityScores []byte
	I             int // index order in the original file
	GCContent     float64
}

// object to hold the input file handles and wrap their close methods
type InputFileReader struct {
	Reader *bufio.Reader
	File   *os.File
	GzFile *os.File
}

func (r *InputFileReader) Close() {
	r.File.Close()
	if r.GzFile != nil {
		r.GzFile.Close()
	}
}

type OutputFileWriter struct {
	File   *os.File
	Writer *gzip.Writer
}

func (w *OutputFileWriter) Close() {
	w.Writer.Flush()
	w.Writer.Close()
	w.File.Close()
}

func GetReader(inputFilepath string) InputFileReader { //(*bufio.Reader, *os.File, *os.File)
	// open the input file for reading
	// the caller needs to run this;
	// defer reader.Close()
	var reader *bufio.Reader
	var file *os.File
	var gzFile *os.File
	bufferSize := 1048576 // default 4096: 4KB ; 1048576 : 1MB ; 10485760 : 10MB

	file, err := os.Open(inputFilepath)
	if err != nil {
		log.Fatalf("Error opening file: %v\n", err)
	}

	if strings.HasSuffix(inputFilepath, ".gz") {
		log.Printf("Opening .gz file %v\n", inputFilepath)
		gz, err := gzip.NewReader(file)
		if err != nil {
			log.Fatalf("Error opening file: %v\n", err)
		}

		reader = bufio.NewReaderSize(gz, bufferSize)
	} else {
		log.Printf("Opening non-compressed file: %v\n", inputFilepath)
		reader = bufio.NewReaderSize(file, bufferSize)
	}

	return InputFileReader{reader, file, gzFile}
}

func GetWriter(outputFilepath string) OutputFileWriter { //(*os.File, *gzip.Writer)
	// initialize the output file writer
	// the caller needs to run this;
	// defer writer.Close()
	outputFile, err := os.Create(outputFilepath)
	if err != nil {
		log.Fatalf("Error creating output file: %v\n", err)
	}

	writer, err := gzip.NewWriterLevel(outputFile, gzip.BestCompression)
	if err != nil {
		log.Fatalf("Error creating output writer: %v\n", err)
	}
	return OutputFileWriter{outputFile, writer}
}

func CalcGCContent(sequence *[]byte) float64 {
	// Calculation for GC content of the DNA sequence
	var numGC float64 = 0.0
	var numBases float64 = float64(len(*sequence))
	for _, base := range *sequence {
		if base == 'G' || base == 'C' {
			numGC = numGC + 1.0
		}
	}
	gcContent := numGC / numBases
	return gcContent
}

func CreateFastqRead(firstLine *[]byte, reader InputFileReader, delim *byte, i *int) FastqRead {
	// reads the next three lines from the reader,
	// and combined with the first line,
	// makes a new FastqRead entry

	sequence, err := reader.Reader.ReadBytes(*delim)
	if err != nil {
		log.Fatalf("Error parsing sequence line in fastq read: %v\n", err)
	}
	plus, err := reader.Reader.ReadBytes(*delim)
	if err != nil {
		log.Fatalf("Error parsing plus line in fastq read: %v\n", err)
	}
	qualityScores, err := reader.Reader.ReadBytes(*delim)
	if err != nil {
		log.Fatalf("Error parsing qualityScores line in fastq read: %v\n", err)
	}

	read := FastqRead{
		Id:            *firstLine,
		Sequence:      sequence,
		Plus:          plus,
		QualityScores: qualityScores,
		I:             *i,
		GCContent:     CalcGCContent(&sequence),
	}
	return read
}

func SortReads(reads *[]FastqRead) {
	// put the read sorting logic in here!
	//
	// sort in-place based on the string of the sequence
	sort.Slice((*reads), func(i, j int) bool { return string((*reads)[i].Sequence) < string((*reads)[j].Sequence) })
}

func SortReadsGC(reads *[]FastqRead) {
	sort.Slice((*reads), func(i, j int) bool { return (*reads)[i].GCContent < (*reads)[j].GCContent })
}

func SortReadsQual(reads *[]FastqRead) {
	sort.Slice((*reads), func(i, j int) bool { return string((*reads)[i].QualityScores) < string((*reads)[j].QualityScores) })
}

func WriteReads(reads *[]FastqRead, writer OutputFileWriter) {
	var n int = 0
	for _, read := range *reads {
		writer.Writer.Write(read.Id)
		writer.Writer.Write(read.Sequence)
		writer.Writer.Write(read.Plus)
		writer.Writer.Write(read.QualityScores)
		n = n + 1
	}
	log.Printf("Wrote %v reads\n", n)
}

func LoadReads(readsBuffer *[]FastqRead, reader InputFileReader, delim *byte) {
	var i int = 0
	for {
		// get the next line
		line, err := reader.Reader.ReadBytes(*delim) // includes the delim in the output line !!
		if err != nil {
			break // end of file io.EOF
		}
		// check if its a FASTQ header line
		if line[0] == '@' {
			i = i + 1
			read := CreateFastqRead(&line, reader, delim, &i)
			*readsBuffer = append(*readsBuffer, read)
		}
	}
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

func SaveOrder(readsBuffer *[]FastqRead) {
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
	reader := GetReader(config.InputFilepath)
	defer reader.Close()

	// output
	writer := GetWriter(config.OutputFilepath)
	defer writer.Close()

	// hold all the reads from the file in here
	reads := []FastqRead{}

	// load all reads from file
	LoadReads(&reads, reader, &config.RecordDelim)
	log.Printf("%v reads loaded\n", len(reads))

	// sort the fastq reads
	SortReads(&reads)
	log.Printf("%v reads after sorting\n", len(reads))

	// write the fastq reads
	log.Printf("Writing to output file %v\n", config.OutputFilepath)
	WriteReads(&reads, writer)

	// save the order of the sorted reads to file
	SaveOrder(&reads)
}

func RunGCSort(config Config) {
	// input
	reader := GetReader(config.InputFilepath)
	defer reader.Close()

	// output
	writer := GetWriter(config.OutputFilepath)
	defer writer.Close()

	// hold all the reads from the file in here
	reads := []FastqRead{}

	// load all reads from file
	LoadReads(&reads, reader, &config.RecordDelim)
	log.Printf("%v reads loaded\n", len(reads))

	// sort the fastq reads
	SortReadsGC(&reads)
	log.Printf("%v reads after sorting\n", len(reads))

	// write the fastq reads
	log.Printf("Writing to output file %v\n", config.OutputFilepath)
	WriteReads(&reads, writer)

	// save the order of the sorted reads to file
	SaveOrder(&reads)
}

func RunQualSort(config Config) {
	// input
	reader := GetReader(config.InputFilepath)
	defer reader.Close()

	// output
	writer := GetWriter(config.OutputFilepath)
	defer writer.Close()

	// hold all the reads from the file in here
	reads := []FastqRead{}

	// load all reads from file
	LoadReads(&reads, reader, &config.RecordDelim)
	log.Printf("%v reads loaded\n", len(reads))

	// sort the fastq reads
	SortReadsQual(&reads)
	log.Printf("%v reads after sorting\n", len(reads))

	// write the fastq reads
	log.Printf("Writing to output file %v\n", config.OutputFilepath)
	WriteReads(&reads, writer)

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
