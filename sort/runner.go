package sort

// wrapper methods to allow for custom handling for different sort methods
// because some sort methods have common setup requirements but others dont

import (
	"code.cloudfoundry.org/bytefmt"
	"log"
	fastq "squish/fastq"
	_io "squish/io"
)

// sort runner method that loads all the fastq reads before sorting
func RunSortLoader(config Config) {
	// run the chosen sorting method on the fastq file

	// input
	reader := _io.GetReader(config.InputFilepath)
	defer reader.Close()

	// output
	writer := _io.GetWriter(config.OutputFilepath)
	defer writer.Close()

	// hold all the reads from the file in here
	reads := []fastq.FastqRead{}

	// load all reads from file
	totalByteSize := fastq.LoadReads(&reads, reader, &config.RecordDelim)
	log.Printf("%v reads loaded (%v)\n", len(reads), bytefmt.ByteSize(uint64(totalByteSize)))

	// sort the fastq reads
	log.Printf("starting read sort")
	config.SortMethod.Func(&reads)
	log.Printf("%v reads after sorting\n", len(reads))

	// write the fastq reads
	log.Printf("Writing to output file %v\n", config.OutputFilepath)
	fastq.WriteReads(&reads, writer)

	// save the order of the sorted reads to file
	fastq.SaveOrder(&reads, config.OrderFilename)
}