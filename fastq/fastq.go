package fastq

// package for handling fastq data

import (
	"bufio"
	"log"
	"os"
	_io "squish/io"
	"strconv"
)

// NOTE: fastq read methods https://pkg.go.dev/github.com/biogo/biogo/io/seqio/fastq#Reader.Read
// https://en.wikipedia.org/wiki/FASTQ_format
type FastqRead struct {
	Id            []byte
	Sequence      []byte
	Plus          []byte
	QualityScores []byte
	I             int // index order in the original file
	GCContent     float64
	Size          int // size in bytes of the read
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

func CreateFastqRead(firstLine *[]byte, reader _io.InputFileReader, delim *byte, i *int) FastqRead {
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
		Size:          len(*firstLine) + len(sequence) + len(plus) + len(qualityScores),
	}
	return read
}

func WriteReads(reads *[]FastqRead, writer _io.OutputFileWriter) {
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

func LoadReads(readsBuffer *[]FastqRead, reader _io.InputFileReader, delim *byte) int {
	var totalSize int = 0
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
			totalSize = totalSize + read.Size
			*readsBuffer = append(*readsBuffer, read)
		}
	}
	return totalSize
}

func SaveOrder(readsBuffer *[]FastqRead, orderFilename string) {
	log.Printf("Saving read order to file %v for %v reads\n", orderFilename, len(*readsBuffer))
	outputFile, err := os.Create(orderFilename)
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
