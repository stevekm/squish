package fastqio

// package for holding the file IO types and methods

import (
	"os"
	// "compress/gzip"
	"bufio"
	"fmt"
	gzip "github.com/klauspost/pgzip"
	"log"
	"strings"
)

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
	reader, err := OpenReader(inputFilepath)
	if err != nil {
		log.Fatalf("Error opening file: %v\n", err)
	}
	return reader
}

func OpenReader(inputFilepath string) (InputFileReader, error) {
	// GetReader hides whether the input is plain FASTQ or gzip-compressed.
	// Callers always receive a buffered reader and only need to defer Close().
	var reader *bufio.Reader
	var file *os.File
	var gzFile *os.File
	bufferSize := 1048576 // default 4096: 4KB ; 1048576 : 1MB ; 10485760 : 10MB

	file, err := os.Open(inputFilepath)
	if err != nil {
		return InputFileReader{}, fmt.Errorf("open input file: %w", err)
	}

	if strings.HasSuffix(inputFilepath, ".gz") {
		log.Printf("Opening .gz file %v\n", inputFilepath)
		gz, err := gzip.NewReader(file)
		if err != nil {
			file.Close()
			return InputFileReader{}, fmt.Errorf("open gzip reader: %w", err)
		}

		reader = bufio.NewReaderSize(gz, bufferSize)
	} else {
		log.Printf("Opening non-compressed file: %v\n", inputFilepath)
		reader = bufio.NewReaderSize(file, bufferSize)
	}

	return InputFileReader{reader, file, gzFile}, nil
}

func GetWriter(outputFilepath string) OutputFileWriter { //(*os.File, *gzip.Writer)
	writer, err := OpenWriter(outputFilepath)
	if err != nil {
		log.Fatalf("Error creating output writer: %v\n", err)
	}
	return writer
}

func OpenWriter(outputFilepath string) (OutputFileWriter, error) {
	// All squish outputs are gzip-compressed FASTQ files. The caller writes raw
	// FASTQ records to Writer and Close flushes the gzip stream and file.
	outputFile, err := os.Create(outputFilepath)
	if err != nil {
		return OutputFileWriter{}, fmt.Errorf("create output file: %w", err)
	}

	writer, err := gzip.NewWriterLevel(outputFile, gzip.BestCompression)
	if err != nil {
		outputFile.Close()
		return OutputFileWriter{}, fmt.Errorf("create gzip writer: %w", err)
	}
	return OutputFileWriter{outputFile, writer}, nil
}
