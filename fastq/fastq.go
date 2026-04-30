package fastq

// package for handling fastq data

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	_io "squish/fastqio"
	"strconv"
	"strings"
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

type ReorderStats struct {
	Reads int
	Bytes int
}

type recordIndexEntry struct {
	Offset int64
	Size   int
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
	read, err := CreateFastqReadE(idOffset, idSize, reader, delim, i, arena)
	if err != nil {
		log.Fatalf("Error creating fastq read: %v\n", err)
	}
	return read
}

func CreateFastqReadE(idOffset int, idSize int, reader _io.InputFileReader, delim *byte, i *int, arena *FastqArena) (FastqRead, error) {
	// The header line has already been read. Pull the remaining three FASTQ
	// lines into the same arena so Record() can later return the exact original
	// four-line record for output.
	sequenceOffset, sequenceSize, err := ReadLineIntoArena(reader, *delim, arena)
	if err != nil {
		return FastqRead{}, fmt.Errorf("parse sequence line in fastq read: %w", err)
	}
	plusOffset, plusSize, err := ReadLineIntoArena(reader, *delim, arena)
	if err != nil {
		return FastqRead{}, fmt.Errorf("parse plus line in fastq read: %w", err)
	}
	qualityScoresOffset, qualityScoresSize, err := ReadLineIntoArena(reader, *delim, arena)
	if err != nil {
		return FastqRead{}, fmt.Errorf("parse qualityScores line in fastq read: %w", err)
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
	return read, nil
}

// ReadNextRead streams one FASTQ record from reader into a small per-record
// arena. The external bucket engine uses this to classify and write records to
// disk without retaining the rest of the input in memory.
func ReadNextRead(reader _io.InputFileReader, delim *byte, i *int) (FastqRead, int, error) {
	return ReadNextReadE(reader, delim, i)
}

func ReadNextReadE(reader _io.InputFileReader, delim *byte, i *int) (FastqRead, int, error) {
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
		read, err := CreateFastqReadE(idOffset, idSize, reader, delim, i, arena)
		if err != nil {
			return FastqRead{}, 0, err
		}
		return read, read.RecordSize, nil
	}
}

func WriteReads(reads *[]FastqRead, writer _io.OutputFileWriter) {
	if err := WriteReadsE(reads, writer); err != nil {
		log.Fatalf("Error writing reads: %v\n", err)
	}
}

func WriteReadsE(reads *[]FastqRead, writer _io.OutputFileWriter) error {
	var n int = 0
	for _, read := range *reads {
		if _, err := writer.Writer.Write(read.Record()); err != nil {
			return fmt.Errorf("write read record: %w", err)
		}
		n = n + 1
	}
	slog.Debug("wrote reads", "count", n)
	return nil
}

func LoadReads(readsBuffer *[]FastqRead, reader _io.InputFileReader, delim *byte) int {
	totalSize, err := LoadReadsE(readsBuffer, reader, delim)
	if err != nil {
		log.Fatalf("Error loading reads: %v\n", err)
	}
	return totalSize
}

func LoadReadsE(readsBuffer *[]FastqRead, reader _io.InputFileReader, delim *byte) (int, error) {
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
			return 0, fmt.Errorf("read fastq header line: %w", err)
		}
		// check if its a FASTQ header line
		line := arena.Data[idOffset : idOffset+idSize]
		if len(line) > 0 && line[0] == '@' {
			i = i + 1
			read, err := CreateFastqReadE(idOffset, idSize, reader, delim, &i, arena)
			if err != nil {
				return 0, err
			}
			totalSize = totalSize + read.RecordSize
			*readsBuffer = append(*readsBuffer, read)
		}
	}
	return totalSize, nil
}

func SaveOrder(readsBuffer *[]FastqRead, orderFilename string) {
	if err := SaveOrderE(readsBuffer, orderFilename); err != nil {
		log.Fatalf("Error saving order: %v\n", err)
	}
}

func SaveOrderE(readsBuffer *[]FastqRead, orderFilename string) error {
	slog.Debug("saving read order", "path", orderFilename, "count", len(*readsBuffer))
	outputFile, err := os.Create(orderFilename)
	if err != nil {
		return fmt.Errorf("create order file: %w", err)
	}
	defer outputFile.Close()

	writer := bufio.NewWriter(outputFile)
	// defer writer.Close()
	defer writer.Flush()
	for _, read := range *readsBuffer {
		_, err := writer.WriteString(strconv.Itoa(read.I) + "\n")
		if err != nil {
			return fmt.Errorf("write order row: %w", err)
		}
	}
	return nil
}

func LoadOrder(orderFilename string) ([]int, error) {
	orderFile, err := os.Open(orderFilename)
	if err != nil {
		return nil, fmt.Errorf("open order file: %w", err)
	}
	defer orderFile.Close()

	order := []int{}
	scanner := bufio.NewScanner(orderFile)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		readIndex, err := strconv.Atoi(line)
		if err != nil {
			return nil, fmt.Errorf("parse order value %q: %w", line, err)
		}
		order = append(order, readIndex)
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scan order file: %w", err)
	}
	return order, nil
}

// NormalizedReadID trims common mate-specific suffixes so R1/R2 headers can be
// compared. It keeps only the first whitespace-delimited field and strips a
// trailing /1 or /2, which covers the common FASTQ header forms.
func NormalizedReadID(id []byte) string {
	value := strings.TrimSpace(string(id))
	value = strings.TrimPrefix(value, "@")
	if fields := strings.Fields(value); len(fields) > 0 {
		value = fields[0]
	}
	if strings.HasSuffix(value, "/1") || strings.HasSuffix(value, "/2") {
		value = value[:len(value)-2]
	}
	return value
}

func LoadNormalizedReadNames(inputFilepath string, delim byte) ([]string, error) {
	reader, err := _io.OpenReader(inputFilepath)
	if err != nil {
		return nil, err
	}
	defer reader.Close()

	names := []string{}
	readIndex := 0
	for {
		read, _, err := ReadNextReadE(reader, &delim, &readIndex)
		if err != nil {
			if errors.Is(err, io.EOF) {
				break
			}
			return nil, fmt.Errorf("read names from fastq: %w", err)
		}
		names = append(names, NormalizedReadID(read.Id()))
	}
	return names, nil
}

// ReorderReadsByOrder applies a sorted R1 order file to a companion FASTQ.
//
// The order file uses one-based original read indexes. For each index in that
// file, the corresponding companion record is emitted to outputFilepath. This
// preserves pairing when R1 defines the sort order.
func ReorderReadsByOrder(
	inputFilepath string,
	outputFilepath string,
	orderFilename string,
	delim byte,
	expectedReads int,
	referenceNames []string,
) (ReorderStats, error) {
	reader, err := _io.OpenReader(inputFilepath)
	if err != nil {
		return ReorderStats{}, err
	}
	defer reader.Close()

	tempDir, err := os.MkdirTemp(filepath.Dir(outputFilepath), ".reorder-*")
	if err != nil {
		return ReorderStats{}, fmt.Errorf("create reorder temp dir: %w", err)
	}
	defer os.RemoveAll(tempDir)

	tempRecordsPath := filepath.Join(tempDir, "records.fastq")
	tempRecords, err := os.Create(tempRecordsPath)
	if err != nil {
		return ReorderStats{}, fmt.Errorf("create reorder temp records: %w", err)
	}

	index := []recordIndexEntry{}
	if expectedReads > 0 {
		index = make([]recordIndexEntry, 0, expectedReads)
	}
	totalBytes := 0
	readIndex := 0
	for {
		read, readSize, err := ReadNextReadE(reader, &delim, &readIndex)
		if err != nil {
			if errors.Is(err, io.EOF) {
				break
			}
			tempRecords.Close()
			return ReorderStats{}, fmt.Errorf("read companion fastq: %w", err)
		}
		if len(referenceNames) > 0 {
			if len(index) >= len(referenceNames) {
				tempRecords.Close()
				return ReorderStats{}, fmt.Errorf("reference read count mismatch for %s: companion has more than %d reads", inputFilepath, len(referenceNames))
			}
			if got, want := NormalizedReadID(read.Id()), referenceNames[len(index)]; got != want {
				tempRecords.Close()
				return ReorderStats{}, fmt.Errorf("read name mismatch for %s at original read %d: got %q, want %q", inputFilepath, len(index)+1, got, want)
			}
		}

		offset, err := tempRecords.Seek(0, io.SeekCurrent)
		if err != nil {
			tempRecords.Close()
			return ReorderStats{}, fmt.Errorf("get temp record offset: %w", err)
		}
		n, err := tempRecords.Write(read.Record())
		if err != nil {
			tempRecords.Close()
			return ReorderStats{}, fmt.Errorf("write temp companion record: %w", err)
		}
		index = append(index, recordIndexEntry{Offset: offset, Size: n})
		totalBytes += readSize
	}
	if err := tempRecords.Close(); err != nil {
		return ReorderStats{}, fmt.Errorf("close temp companion records: %w", err)
	}

	if expectedReads > 0 && len(index) != expectedReads {
		return ReorderStats{}, fmt.Errorf("read count mismatch for %s: got %d, want %d", inputFilepath, len(index), expectedReads)
	}
	if len(referenceNames) > 0 && len(referenceNames) != len(index) {
		return ReorderStats{}, fmt.Errorf("reference read count mismatch for %s: got %d reference names for %d reads", inputFilepath, len(referenceNames), len(index))
	}

	orderFile, err := os.Open(orderFilename)
	if err != nil {
		return ReorderStats{}, fmt.Errorf("open order file: %w", err)
	}
	defer orderFile.Close()

	tempRecordsReader, err := os.Open(tempRecordsPath)
	if err != nil {
		return ReorderStats{}, fmt.Errorf("open temp companion records: %w", err)
	}
	defer tempRecordsReader.Close()

	seen := make([]bool, len(index))
	writer, err := _io.OpenWriter(outputFilepath)
	if err != nil {
		return ReorderStats{}, err
	}
	defer writer.Close()

	scanner := bufio.NewScanner(orderFile)
	outputIndex := 0
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		originalIndex, err := strconv.Atoi(line)
		if err != nil {
			return ReorderStats{}, fmt.Errorf("parse order value %q: %w", line, err)
		}
		outputIndex++
		if originalIndex < 1 || originalIndex > len(index) {
			return ReorderStats{}, fmt.Errorf("order row %d references read %d outside range [1, %d]", outputIndex, originalIndex, len(index))
		}
		if seen[originalIndex-1] {
			return ReorderStats{}, fmt.Errorf("order row %d repeats read %d", outputIndex, originalIndex)
		}
		seen[originalIndex-1] = true
		recordIndex := index[originalIndex-1]
		record := make([]byte, recordIndex.Size)
		if _, err := tempRecordsReader.ReadAt(record, recordIndex.Offset); err != nil {
			return ReorderStats{}, fmt.Errorf("read temp companion record %d: %w", originalIndex, err)
		}
		if _, err := writer.Writer.Write(record); err != nil {
			return ReorderStats{}, fmt.Errorf("write reordered read %d: %w", originalIndex, err)
		}
	}
	if err := scanner.Err(); err != nil {
		return ReorderStats{}, fmt.Errorf("scan order file: %w", err)
	}
	if outputIndex != len(index) {
		return ReorderStats{}, fmt.Errorf("order length mismatch for %s: got %d order rows for %d reads", inputFilepath, outputIndex, len(index))
	}

	return ReorderStats{Reads: len(index), Bytes: totalBytes}, nil
}
