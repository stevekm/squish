package io

// package for holding the file IO types and methods

import (
	"os"
	// "compress/gzip"
	"bufio"
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
