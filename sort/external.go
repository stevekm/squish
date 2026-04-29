package sort

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	go_sort "sort"
	"strconv"

	"code.cloudfoundry.org/bytefmt"
	fastq "squish/fastq"
	_io "squish/io"
)

// ExternalBucketConfig contains the file paths and parsing settings needed by
// the file-backed sorting engine.
type ExternalBucketConfig struct {
	InputFilepath  string
	OutputFilepath string
	OrderFilepath  string
	TempDir        string
	RecordDelim    byte
}

// bucketWriter owns the temporary FASTQ bucket and a sidecar order file.
//
// The order file records each read's original input index because temporary
// FASTQ files only contain the four FASTQ lines. Restoring that index after a
// bucket is reloaded keeps sort tie-breaks and order.txt output correct.
type bucketWriter struct {
	file        *os.File
	writer      *bufio.Writer
	path        string
	orderFile   *os.File
	orderWriter *bufio.Writer
	orderPath   string
}

// maxOpenBucketWriters caps file descriptor usage when a bucket strategy has a
// large number of possible buckets. Closed buckets are reopened in append mode
// if more reads later map to the same bucket.
const maxOpenBucketWriters = 64

// RunExternalBucketSort streams the input into temporary bucket files, then
// sorts and emits one bucket at a time.
//
// Memory usage is bounded by the largest single bucket instead of the full
// input file. Ordered bucket strategies produce exact global sorts by writing
// sorted buckets in bucket ID order.
func RunExternalBucketSort(config ExternalBucketConfig, sorter SortStrategy, bucketer BucketStrategy) error {
	if !bucketer.OrderedFor(sorter) {
		slog.Debug(
			"bucket strategy is not globally ordered for sorter; output will be deterministic but not a strict global sort",
			"sorter", sorter.Name(),
			"bucketer", bucketer.Name(),
		)
	}

	// Start from a clean temp directory so append-mode bucket files cannot pick
	// up stale records from a previous failed or interrupted run.
	if err := os.RemoveAll(config.TempDir); err != nil {
		return fmt.Errorf("clean temp dir: %w", err)
	}
	if err := os.MkdirAll(config.TempDir, 0755); err != nil {
		return fmt.Errorf("create temp dir: %w", err)
	}

	bucketPaths, bucketOrderPaths, bucketSizes, totalReads, totalBytes, err := writeBuckets(config, bucketer)
	if err != nil {
		return err
	}
	slog.Info(
		"external buckets written",
		"reads", totalReads,
		"size", bytefmt.ByteSize(uint64(totalBytes)),
		"buckets_used", len(bucketPaths),
		"bucketer", bucketer.Name(),
	)

	if err := sortBucketsToOutput(config, sorter, bucketer, bucketPaths, bucketOrderPaths, bucketSizes); err != nil {
		return err
	}

	return nil
}

// writeBuckets is the streaming phase. It reads one FASTQ record at a time,
// computes the bucket ID, and appends the raw record to that bucket's temp file.
func writeBuckets(config ExternalBucketConfig, bucketer BucketStrategy) (map[int]string, map[int]string, map[int]int64, int, int, error) {
	reader := _io.GetReader(config.InputFilepath)
	defer reader.Close()

	writers := map[int]*bucketWriter{}
	bucketPaths := map[int]string{}
	bucketOrderPaths := map[int]string{}
	bucketSizes := map[int]int64{}
	totalReads := 0
	totalBytes := 0
	readIndex := 0

	defer func() {
		for _, bucket := range writers {
			closeBucketWriter(bucket)
		}
	}()

	for {
		// ReadNextRead returns a small one-record arena, so this loop does not
		// retain the full input in memory during bucketing.
		read, readSize, err := fastq.ReadNextRead(reader, &config.RecordDelim, &readIndex)
		if err != nil {
			if errors.Is(err, io.EOF) {
				break
			}
			return nil, nil, nil, 0, 0, fmt.Errorf("read fastq record: %w", err)
		}

		bucketID := bucketer.BucketID(read)
		if bucketID < 0 || bucketID >= bucketer.BucketCount() {
			return nil, nil, nil, 0, 0, fmt.Errorf("bucket id %d outside range [0, %d)", bucketID, bucketer.BucketCount())
		}

		bucket, ok := writers[bucketID]
		if !ok {
			if len(writers) >= maxOpenBucketWriters {
				// Close any currently open bucket to stay under the descriptor
				// cap. The bucket will be reopened in append mode when needed.
				closeOneBucketWriter(writers)
			}
			var err error
			bucket, err = openBucketWriter(config.TempDir, bucketID)
			if err != nil {
				return nil, nil, nil, 0, 0, err
			}
			writers[bucketID] = bucket
			bucketPaths[bucketID] = bucket.path
			bucketOrderPaths[bucketID] = bucket.orderPath
		}

		n, err := bucket.writer.Write(read.Record())
		if err != nil {
			return nil, nil, nil, 0, 0, fmt.Errorf("write bucket %d: %w", bucketID, err)
		}
		if _, err := bucket.orderWriter.WriteString(strconv.Itoa(read.I) + "\n"); err != nil {
			return nil, nil, nil, 0, 0, fmt.Errorf("write bucket order %d: %w", bucketID, err)
		}
		bucketSizes[bucketID] += int64(n)
		totalReads++
		totalBytes += readSize
	}

	return bucketPaths, bucketOrderPaths, bucketSizes, totalReads, totalBytes, nil
}

// openBucketWriter opens both the FASTQ bucket and its order sidecar. Files are
// opened with O_APPEND so buckets can be closed and reopened safely.
func openBucketWriter(tempDir string, bucketID int) (*bucketWriter, error) {
	path := filepath.Join(tempDir, fmt.Sprintf("bucket-%06d.fastq", bucketID))
	file, err := os.OpenFile(path, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0644)
	if err != nil {
		return nil, fmt.Errorf("create bucket %d: %w", bucketID, err)
	}
	orderPath := filepath.Join(tempDir, fmt.Sprintf("bucket-%06d.order", bucketID))
	orderFile, err := os.OpenFile(orderPath, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0644)
	if err != nil {
		file.Close()
		return nil, fmt.Errorf("create bucket order %d: %w", bucketID, err)
	}
	return &bucketWriter{
		file:        file,
		writer:      bufio.NewWriter(file),
		path:        path,
		orderFile:   orderFile,
		orderWriter: bufio.NewWriter(orderFile),
		orderPath:   orderPath,
	}, nil
}

// closeOneBucketWriter evicts one open bucket writer from the open-writer map.
// Any bucket is acceptable because future writes reopen the file in append mode.
func closeOneBucketWriter(writers map[int]*bucketWriter) {
	for bucketID, bucket := range writers {
		closeBucketWriter(bucket)
		delete(writers, bucketID)
		return
	}
}

// closeBucketWriter flushes both buffered files before closing them. Errors are
// logged because this helper is also used from cleanup paths.
func closeBucketWriter(bucket *bucketWriter) {
	if err := bucket.writer.Flush(); err != nil {
		slog.Debug("could not flush bucket", "path", bucket.path, "error", err)
	}
	if err := bucket.orderWriter.Flush(); err != nil {
		slog.Debug("could not flush bucket order", "path", bucket.orderPath, "error", err)
	}
	if err := bucket.file.Close(); err != nil {
		slog.Debug("could not close bucket", "path", bucket.path, "error", err)
	}
	if err := bucket.orderFile.Close(); err != nil {
		slog.Debug("could not close bucket order", "path", bucket.orderPath, "error", err)
	}
}

// sortBucketsToOutput is the bounded-memory sort phase. Each bucket is loaded
// into the arena representation, sorted in memory, appended to the gzip output,
// and then deleted.
func sortBucketsToOutput(
	config ExternalBucketConfig,
	sorter SortStrategy,
	bucketer BucketStrategy,
	bucketPaths map[int]string,
	bucketOrderPaths map[int]string,
	bucketSizes map[int]int64,
) error {
	outputWriter := _io.GetWriter(config.OutputFilepath)
	defer outputWriter.Close()

	orderFile, err := os.Create(config.OrderFilepath)
	if err != nil {
		return fmt.Errorf("create order file: %w", err)
	}
	defer orderFile.Close()
	orderWriter := bufio.NewWriter(orderFile)
	defer orderWriter.Flush()

	for bucketID := 0; bucketID < bucketer.BucketCount(); bucketID++ {
		bucketPath, ok := bucketPaths[bucketID]
		if !ok {
			continue
		}

		reads, err := loadBucket(bucketPath, bucketOrderPaths[bucketID], config.RecordDelim)
		if err != nil {
			return err
		}
		go_sort.Slice(reads, func(i, j int) bool {
			return sorter.Less(reads[i], reads[j])
		})

		// Ordered bucket strategies rely on this append order: bucket 0 first,
		// then bucket 1, and so on. Each bucket is already internally sorted.
		for _, read := range reads {
			if _, err := outputWriter.Writer.Write(read.Record()); err != nil {
				return fmt.Errorf("write output record: %w", err)
			}
			if _, err := orderWriter.WriteString(strconv.Itoa(read.I) + "\n"); err != nil {
				return fmt.Errorf("write order record: %w", err)
			}
		}

		slog.Debug(
			"external bucket sorted",
			"bucket", bucketID,
			"reads", len(reads),
			"size", bytefmt.ByteSize(uint64(bucketSizes[bucketID])),
		)

		if err := os.Remove(bucketPath); err != nil {
			return fmt.Errorf("remove bucket %d: %w", bucketID, err)
		}
		if err := os.Remove(bucketOrderPaths[bucketID]); err != nil {
			return fmt.Errorf("remove bucket order %d: %w", bucketID, err)
		}
	}

	return nil
}

// loadBucket reloads one temporary FASTQ bucket into memory and restores the
// original global read indexes from the sidecar order file.
func loadBucket(bucketPath string, orderPath string, delim byte) ([]fastq.FastqRead, error) {
	reader := _io.GetReader(bucketPath)
	defer reader.Close()

	reads := []fastq.FastqRead{}
	fastq.LoadReads(&reads, reader, &delim)
	if err := restoreOriginalOrder(reads, orderPath); err != nil {
		return nil, err
	}
	return reads, nil
}

// restoreOriginalOrder copies the global input order back onto reads after a
// temp bucket has been loaded. LoadReads assigns per-bucket indexes, so this
// sidecar is required for stable tie-breaks and correct final order.txt output.
func restoreOriginalOrder(reads []fastq.FastqRead, orderPath string) error {
	orderFile, err := os.Open(orderPath)
	if err != nil {
		return fmt.Errorf("open bucket order file: %w", err)
	}
	defer orderFile.Close()

	scanner := bufio.NewScanner(orderFile)
	i := 0
	for scanner.Scan() {
		if i >= len(reads) {
			return fmt.Errorf("bucket order file has more rows than reads")
		}
		order, err := strconv.Atoi(scanner.Text())
		if err != nil {
			return fmt.Errorf("parse bucket order: %w", err)
		}
		reads[i].I = order
		i++
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("scan bucket order: %w", err)
	}
	if i != len(reads) {
		return fmt.Errorf("bucket order file has %d rows for %d reads", i, len(reads))
	}
	return nil
}

// MustRunExternalBucketSort is a CLI-oriented wrapper around
// RunExternalBucketSort that exits on error.
func MustRunExternalBucketSort(config ExternalBucketConfig, sorter SortStrategy, bucketer BucketStrategy) {
	if err := RunExternalBucketSort(config, sorter, bucketer); err != nil {
		log.Fatalf("ERROR: external bucket sort failed: %v\n", err)
	}
}
