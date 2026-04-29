package fastq

// package for handling fastq data

import (
	"bufio"
	"errors"
	"io"
	"log"
	"log/slog"
	"os"
	_io "squish/io"
	"strconv"
)

// NOTE: fastq read methods https://pkg.go.dev/github.com/biogo/biogo/io/seqio/fastq#Reader.Read
// https://en.wikipedia.org/wiki/FASTQ_format
type FastqRead struct {
	Arena              *FastqArena
	RecordOffset       int
	RecordSize         int
	IdOffset           int
	IdSize             int
	SequenceOffset     int
	SequenceSize       int
	PlusOffset         int
	PlusSize           int
	QualityScoreOffset int
	QualityScoreSize   int
	I                  int // index order in the original file
	GCContent          float64
}

type FastqArena struct {
	Data []byte
}

func (a *FastqArena) Append(line []byte) (int, int) {
	offset := len(a.Data)
	a.Data = append(a.Data, line...)
	return offset, len(line)
}

func (read FastqRead) Id() []byte {
	return read.Arena.Data[read.IdOffset : read.IdOffset+read.IdSize]
}

func (read FastqRead) Sequence() []byte {
	return read.Arena.Data[read.SequenceOffset : read.SequenceOffset+read.SequenceSize]
}

func (read FastqRead) Plus() []byte {
	return read.Arena.Data[read.PlusOffset : read.PlusOffset+read.PlusSize]
}

func (read FastqRead) QualityScores() []byte {
	return read.Arena.Data[read.QualityScoreOffset : read.QualityScoreOffset+read.QualityScoreSize]
}

func (read FastqRead) Record() []byte {
	return read.Arena.Data[read.RecordOffset : read.RecordOffset+read.RecordSize]
}

func CalcGCContent(sequence []byte) float64 {
	// Calculation for GC content of the DNA sequence
	var numGC float64 = 0.0
	var numBases float64 = float64(len(sequence))
	for _, base := range sequence {
		if base == 'G' || base == 'C' {
			numGC = numGC + 1.0
		}
	}
	if numBases == 0 {
		return 0.0
	}
	gcContent := numGC / numBases
	return gcContent
}

func ReadLineIntoArena(reader _io.InputFileReader, delim byte, arena *FastqArena) (int, int, error) {
	offset := len(arena.Data)
	for {
		line, err := reader.Reader.ReadSlice(delim)
		if len(line) > 0 {
			arena.Data = append(arena.Data, line...)
		}
		if err == nil {
			return offset, len(arena.Data) - offset, nil
		}
		if errors.Is(err, bufio.ErrBufferFull) {
			continue
		}
		if len(arena.Data) > offset {
			if errors.Is(err, io.EOF) {
				return offset, len(arena.Data) - offset, nil
			}
			return offset, len(arena.Data) - offset, err
		}
		return 0, 0, err
	}
}

func CreateFastqRead(idOffset int, idSize int, reader _io.InputFileReader, delim *byte, i *int, arena *FastqArena) FastqRead {
	// reads the next three lines from the reader,
	// and combined with the first line,
	// makes a new FastqRead entry
	sequenceOffset, sequenceSize, err := ReadLineIntoArena(reader, *delim, arena)
	if err != nil {
		log.Fatalf("Error parsing sequence line in fastq read: %v\n", err)
	}
	plusOffset, plusSize, err := ReadLineIntoArena(reader, *delim, arena)
	if err != nil {
		log.Fatalf("Error parsing plus line in fastq read: %v\n", err)
	}
	qualityScoresOffset, qualityScoresSize, err := ReadLineIntoArena(reader, *delim, arena)
	if err != nil {
		log.Fatalf("Error parsing qualityScores line in fastq read: %v\n", err)
	}

	recordOffset := idOffset
	recordSize := idSize + sequenceSize + plusSize + qualityScoresSize
	read := FastqRead{
		Arena:              arena,
		RecordOffset:       recordOffset,
		RecordSize:         recordSize,
		IdOffset:           idOffset,
		IdSize:             idSize,
		SequenceOffset:     sequenceOffset,
		SequenceSize:       sequenceSize,
		PlusOffset:         plusOffset,
		PlusSize:           plusSize,
		QualityScoreOffset: qualityScoresOffset,
		QualityScoreSize:   qualityScoresSize,
		I:                  *i,
		GCContent:          CalcGCContent(arena.Data[sequenceOffset : sequenceOffset+sequenceSize]),
	}
	return read
}

func WriteReads(reads *[]FastqRead, writer _io.OutputFileWriter) {
	var n int = 0
	for _, read := range *reads {
		writer.Writer.Write(read.Record())
		n = n + 1
	}
	slog.Debug("wrote reads", "count", n)
}

func LoadReads(readsBuffer *[]FastqRead, reader _io.InputFileReader, delim *byte) int {
	var totalSize int = 0
	var i int = 0
	arena := &FastqArena{}
	for {
		// get the next line
		idOffset, idSize, err := ReadLineIntoArena(reader, *delim, arena)
		if err != nil {
			if errors.Is(err, io.EOF) {
				break
			}
			log.Fatalf("Error reading fastq header line: %v\n", err)
		}
		// check if its a FASTQ header line
		line := arena.Data[idOffset : idOffset+idSize]
		if len(line) > 0 && line[0] == '@' {
			i = i + 1
			read := CreateFastqRead(idOffset, idSize, reader, delim, &i, arena)
			totalSize = totalSize + read.RecordSize
			*readsBuffer = append(*readsBuffer, read)
		}
	}
	return totalSize
}

func SaveOrder(readsBuffer *[]FastqRead, orderFilename string) {
	slog.Debug("saving read order", "path", orderFilename, "count", len(*readsBuffer))
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
