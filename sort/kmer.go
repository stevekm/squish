package sort

// parsing for KMers for KMer sorting

import (
	"fmt"
	"log"
	"io"
	fastq "squish/fastq"
	"code.cloudfoundry.org/bytefmt"
)

// Pair represents a k-mer and its associated read
type Pair struct {
	Kmer string
	Read *fastq.FastqRead
}

func SortKmer(reads *[]fastq.FastqRead, config Config) {
	log.Println("---- USING KMER SORT ----")
	// reads = nil

	// input
	fastqReader := fastq.NewFastqReader(config.InputFilepath)
	defer fastqReader.Reader.Close()

	// reads := []fastq.FastqRead{}

	for {
		read, err := fastqReader.Next()
		if err == io.EOF {
			break

		}
		if err != nil {
			// TODO: do something about error here
			log.Fatalf("Error reading lines from fastq file: %v", err)
		}

		*reads = append(*reads, read)
	}
	log.Printf("%v reads loaded (%v)\n", fastqReader.NRead, bytefmt.ByteSize(uint64(fastqReader.BytesRead)))

	// // output
	// writer := _io.GetWriter(config.OutputFilepath)
	// defer writer.Close()

	kmersList := []Pair{}

	for _, read := range *reads {
		kmers := ExtractKmers(string(read.Sequence))
		// fmt.Printf("READ: %v, KMERS: %v\n", string(read.Sequence), kmers)
		for _, kmer := range kmers {
			kmersList = append(kmersList, Pair{kmer, &read})
		}
	}

	sortedKmersList := RadixSortKmers(kmersList)

	for _, kmer := range sortedKmersList {
		fmt.Printf("KMER: %v, READ: %v\n", kmer.Kmer, string(kmer.Read.Sequence))
	}


	log.Println("---- DONE KMER SORT ----")
}

func ExtractKmers(read string) []string {
	k := 3 // Assuming k-mer length of 3
	var kmers []string
	for i := 0; i <= len(read)-k; i++ {
		kmers = append(kmers, read[i:i+k])
	}
	return kmers
}

// RadixSortKmers sorts k-mers using radix sort
func RadixSortKmers(kmerReadPairs []Pair) []Pair {
	// Determine the maximum length of k-mers
	maxLength := 0
	for _, pair := range kmerReadPairs {
		if len(pair.Kmer) > maxLength {
			maxLength = len(pair.Kmer)
		}
	}

	// Perform radix sort for each character position
	for i := maxLength - 1; i >= 0; i-- {
		// Use counting sort to distribute k-mers into buckets based on the i-th character
		buckets := make([][]Pair, 256) // Assuming ASCII characters
		for _, pair := range kmerReadPairs {
			var index int
			if i < len(pair.Kmer) {
				index = int(pair.Kmer[i])
			}
			buckets[index] = append(buckets[index], pair)
		}

		// Combine the buckets to form a new order of k-mers and associated reads
		// sortMethodsDescr := map[string]string{}
		// kmerReadPairsMap := map[*fastq.FastqRead]bool{}
		kmerReadPairs = nil

		for _, bucket := range buckets {
			// TODO: this is outputting all the reads multiple times !! Fix this

			// $ ./squish -m kmer data/test1.fastq.gz test1.squish.fastq.gz | wc -l
			// 2023/11/30 09:44:16 3276 reads loaded (854.2K)
			// 658476

			kmerReadPairs = append(kmerReadPairs, bucket...)

			// only output unique reads
			// TODO this does not seem to work right it only outputs a few reads
			// for _, pair := range bucket {
			// 	// check if the read has already been processed
			// 	_, ok := kmerReadPairsMap[pair.Read]
			// 	if !ok {
			// 		kmerReadPairs = append(kmerReadPairs, pair)
			// 		kmerReadPairsMap[pair.Read] = true
			// 	}
			// }

		}
	}

	return kmerReadPairs
}
