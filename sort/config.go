package sort

// put the configuration methods for sorting here

import (
	"fmt"
	"log"
	fastq "squish/fastq"
	"time"
)

// default values; make sure these match!
const DefaultSortMethod string = "alpha"
const DefaultSortDescrption string = "Alphabetical sort on sequence"

// holds all the parameters that are used for sorting
type Config struct {
	SortMethod       SortMethod
	InputFilepath    string
	InputFileSize    int64
	OutputFilepath   string
	// RecordDelim      byte
	// RecordHeaderChar byte
	TimeStart        time.Time
	OrderFilename    string
}

func (config *Config) Run() {
	log.Printf("Using sort method: %v\n", config.SortMethod.CLIArg)
	config.SortMethod.Runner(*config)
}

// keep the CLIarg's associated with sort method and runner methods
type SortMethod struct {
	CLIArg      string
	Description string
	Func        func(*[]fastq.FastqRead, Config)
	Runner      func(Config)
}

// parse & return the available sorting methods
func GetSortingMethods() (map[string]SortMethod, string) {
	sortMethodMap := map[string]SortMethod{
		// alpha sort by default!
		DefaultSortMethod: SortMethod{
			CLIArg:      DefaultSortMethod,
			Description: DefaultSortDescrption,
			Func:        SortReadsSequence,
			Runner:      RunSortLoader}, // NOTE: MAKE SURE THE DEFAULT LOADER IS CORRECT FOR THE DEFAULT METHOD!
		"gc": SortMethod{
			CLIArg:      "gc",
			Description: "GC Content Sort",
			Func:        SortReadsGC,
			Runner:      RunSortLoader},
		"qual": SortMethod{
			CLIArg:      "qual",
			Description: "Quality score sort",
			Func:        SortReadsQual,
			Runner:      RunSortLoader},
		"alpha-heap": SortMethod{
			CLIArg:      "alpha-heap",
			Description: "Sequence alpha heap sort",
			Func:        HeapSortSequence,
			Runner:      RunSortLoader},
		"kmer": SortMethod{
			CLIArg:      "kmer",
			Description: "Kmer sort",
			Func:        SortKmer,
			Runner:      RunSort},
	}
	// minimal map for help text printing
	sortMethodsDescr := map[string]string{}
	for key, value := range sortMethodMap {
		sortMethodsDescr[key] = value.Description
	}
	// help text
	// TODO: make this cleaner for printing, its kinda hard to read in the console
	sortMethodOptionStr := fmt.Sprintf("Options: %v", sortMethodsDescr)
	return sortMethodMap, sortMethodOptionStr
}
