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

// FastqRead is a lightweight view into a FastqArena.
//
// Instead of allocating four separate byte slices per read, each field stores
// an offset and size into the shared arena buffer. That keeps the read structs
// cheap to copy during sorting while preserving fast access to the original
// FASTQ record bytes.
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

// FastqArena owns the raw FASTQ bytes referenced by one or more FastqRead
// values. In memory mode there is one arena for the full input; in external
// bucket mode each streamed record or loaded bucket gets a smaller arena.
type FastqArena struct {
	Data []byte
}

// Append adds a line or record fragment to the arena and returns the range
// that was appended. Callers store those ranges on FastqRead instead of
// retaining separate slices.
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

// ReadLineIntoArena reads one delimited FASTQ line into the arena.
//
// bufio.Reader.ReadSlice can return ErrBufferFull for long lines. The loop
// keeps appending chunks until it sees the delimiter or EOF, so the parser can
// handle reads longer than the bufio internal buffer.
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
	// The header line has already been read. Pull the remaining three FASTQ
	// lines into the same arena so Record() can later return the exact original
	// four-line record for output.
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

// ReadNextRead streams one FASTQ record from reader into a small per-record
// arena. The external bucket engine uses this to classify and write records to
// disk without retaining the rest of the input in memory.
func ReadNextRead(reader _io.InputFileReader, delim *byte, i *int) (FastqRead, int, error) {
	arena := &FastqArena{}
	for {
		idOffset, idSize, err := ReadLineIntoArena(reader, *delim, arena)
		if err != nil {
			return FastqRead{}, 0, err
		}

		line := arena.Data[idOffset : idOffset+idSize]
		if len(line) == 0 || line[0] != '@' {
			// Skip non-header lines and reset the arena so stray bytes are not
			// retained in the next candidate record.
			arena = &FastqArena{}
			continue
		}

		// CreateFastqRead reads the following sequence, plus, and quality lines
		// into the same per-record arena.
		*i = *i + 1
		read := CreateFastqRead(idOffset, idSize, reader, delim, i, arena)
		return read, read.RecordSize, nil
	}
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
	// LoadReads is the memory-engine parser: all records share one arena and
	// the returned FastqRead structs only carry offsets into that buffer.
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
