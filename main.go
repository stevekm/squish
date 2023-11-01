package main

import (
	// "fmt"
	// "compress/gzip"
	"bufio"
	gzip "github.com/klauspost/pgzip"
	"log"
	"os"
	"sort"
	"runtime/pprof"
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

func CreateFastqRead(firstLine *[]byte, reader *bufio.Reader, delim *byte) FastqRead {
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
	for _, read := range *reads {
		writer.Write(read.Id)
		writer.Write(read.Sequence)
		writer.Write(read.Plus)
		writer.Write(read.QualityScores)
	}
}

func LoadReads(readsBuffer *[]FastqRead, reader *bufio.Reader, delim *byte) {
	for {
		// get the next line
		line, err := reader.ReadBytes(*delim) // includes the delim in the output line !!
		if err != nil {
			break // end of file io.EOF
		}
		// check if its a FASTQ header line
		if line[0] == '@' {
			read := CreateFastqRead(&line, reader, delim)
			*readsBuffer = append(*readsBuffer, read)
		}
	}
}

func Profiling()(*os.File, *os.File){
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

func main() {
	// start profiler
	cpuFile, memFile := Profiling()
	defer cpuFile.Close()
	defer memFile.Close()
	defer pprof.StopCPUProfile()

	// get command line args
	inputFilepath := os.Args[1]
	outputFilepath := os.Args[2]
	var delim byte = '\n'

	inputFileSize, err := GetFileSize(inputFilepath)
	if err != nil {
		log.Printf("WARNING: could not get size for file %v\n", inputFilepath)
	}
	log.Printf("Input file %v of size %v Bytes\n", inputFilepath, inputFileSize)

	// input
	reader, file, gzFile := GetReader(inputFilepath)
	defer file.Close()
	defer gzFile.Close()

	// output
	outputFile, writer := GetWriter(outputFilepath)
	defer outputFile.Close()

	// hold all the reads from the file in here
	reads := []FastqRead{}

	// load all reads from file
	LoadReads(&reads, reader, &delim)
	numReads := len(reads)
	log.Printf("%v reads loaded\n", numReads)

	// sort the fastq reads
	SortReads(&reads)
	log.Printf("%v reads sorted\n", len(reads))

	pprof.WriteHeapProfile(memFile)

	// write the fastq reads
	WriteReads(&reads, writer)
	writer.Flush()
	writer.Close()

	// print some stuff to the console log
	outputFileSize, err := GetFileSize(outputFilepath)
	if err != nil {
		log.Printf("WARNING: could not get size for file %v\n", outputFilepath)
	}
	sizeDifference := inputFileSize - outputFileSize
	log.Printf("Output file created %v of size %v Bytes\n", outputFilepath, outputFileSize)
	log.Printf("Size reduced by %v Bytes (%.4f)\n",
		sizeDifference,
		float64(sizeDifference) / float64(inputFileSize),
	)
}
