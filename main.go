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

func GetReader(inputFilepath string) (*bufio.Reader, *os.File, *os.File) {
	// open the input file for reading
	// the caller needs to run this;
	// defer file.Close()
	// defer gzFile.Close()
	var reader *bufio.Reader
	var file *os.File
	var gzFile *os.File
	bufferSize := 1048576 // default 4096: 4KB ; 1048576 : 1MB ; 10485760 : 10MB

	file, err := os.Open(inputFilepath)
	if err != nil {
		log.Fatalf("Error opening file: %v\n", err)
	}

	gz, err := gzip.NewReader(file)
	if err != nil {
		log.Fatalf("Error opening file: %v\n", err)
	}

	reader = bufio.NewReaderSize(gz, bufferSize)

	return reader, file, gzFile
}

func GetWriter(outputFilepath string) (*os.File, *gzip.Writer) {
	// initialize the output file writer
	// the caller needs to run this;
	// GzWriter.Flush()
	// GzWriter.Close()
	// File.Close()
	outputFile, err := os.Create(outputFilepath)
	if err != nil {
		log.Fatalf("Error creating output file: %v\n", err)
	}

	writer, err := gzip.NewWriterLevel(outputFile, gzip.BestCompression)
	if err != nil {
		log.Fatalf("Error creating output writer: %v\n", err)
	}
	return outputFile, writer
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

func CreateFastqRead(firstLine *[]byte, reader *bufio.Reader, delim *byte, i *int) FastqRead {
	// reads the next three lines from the reader,
	// and combined with the first line,
	// makes a new FastqRead entry

	sequence, err := reader.ReadBytes(*delim)
	if err != nil {
		log.Fatalf("Error parsing sequence line in fastq read: %v\n", err)
	}
	plus, err := reader.ReadBytes(*delim)
	if err != nil {
		log.Fatalf("Error parsing plus line in fastq read: %v\n", err)
	}
	qualityScores, err := reader.ReadBytes(*delim)
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

func WriteReads(reads *[]FastqRead, writer *gzip.Writer) {
	var n int = 0
	for _, read := range *reads {
		writer.Write(read.Id)
		writer.Write(read.Sequence)
		writer.Write(read.Plus)
		writer.Write(read.QualityScores)
		n = n + 1
	}
	log.Printf("Wrote %v reads\n", n)
	writer.Flush()
	writer.Close()
}

func LoadReads(readsBuffer *[]FastqRead, reader *bufio.Reader, delim *byte) {
	var i int = 0
	for {
		// get the next line
		line, err := reader.ReadBytes(*delim) // includes the delim in the output line !!
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

func Profiling() (*os.File, *os.File) {
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
	// pprof.WriteHeapProfile(memFile)
	// defer pprof.StopCPUProfile()
	//
	// see also; https://pkg.go.dev/net/http/pprof

	// Start CPU profiling
	cpuFile, err := os.Create("cpu.prof")
	if err != nil {
		log.Fatal(err)
	}
	pprof.StartCPUProfile(cpuFile)

	// Start memory profiling file
	memFile, err := os.Create("mem.prof")
	if err != nil {
		log.Fatal(err)
	}

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

type Config struct {
	SortMethod       string
	InputFilepath    string
	OutputFilepath   string
	RecordDelim      byte
	RecordHeaderChar byte
	TimeStart        time.Time
}

func main() {
	timeStart := time.Now()
	// start profiler
	cpuFile, memFile := Profiling()
	defer cpuFile.Close()
	defer memFile.Close()
	defer pprof.StopCPUProfile()

	// get command line args
	printVersion := flag.Bool("v", false, "print version information")
	flag.Parse()
	cliArgs := flag.Args() // all positional args passed
	if *printVersion {
		PrintVersionAndQuit()
	}
	MinCliPosArgs(cliArgs, 2)
	inputFilepath := cliArgs[0]
	outputFilepath := cliArgs[1]
	var delim byte = '\n'

	config := Config{
		SortMethod:       "alpha",
		InputFilepath:    inputFilepath,
		OutputFilepath:   outputFilepath,
		RecordDelim:      delim,
		RecordHeaderChar: '@',
		TimeStart:        timeStart,
	}

	//
	//
	//

	inputFileSize := LogFileSize(config.InputFilepath, "Input")

	// input
	reader, file, gzFile := GetReader(config.InputFilepath)
	defer file.Close()
	defer gzFile.Close()

	// output
	outputFile, writer := GetWriter(outputFilepath)
	defer outputFile.Close()

	// hold all the reads from the file in here
	reads := []FastqRead{}

	// load all reads from file
	LoadReads(&reads, reader, &delim)
	log.Printf("%v reads loaded\n", len(reads))

	// sort the fastq reads
	SortReads(&reads)
	log.Printf("%v reads after sorting\n", len(reads))

	// write the fastq reads
	log.Printf("Writing to output file %v\n", outputFilepath)
	WriteReads(&reads, writer)

	pprof.WriteHeapProfile(memFile)

	SaveOrder(&reads)

	//
	//
	//

	// print some stuff to the console log
	timeStop := time.Now()
	timeDuration := timeStop.Sub(timeStart)
	outputFileSize := LogFileSize(config.OutputFilepath, "Output")
	sizeDifference := inputFileSize - outputFileSize
	sizeDifferenceBytes := bytefmt.ByteSize(uint64(sizeDifference))

	log.Printf("Size reduced by %v Bytes (%.4f) in %v\n",
		sizeDifferenceBytes,
		float64(sizeDifference)/float64(inputFileSize),
		timeDuration,
	)
}
