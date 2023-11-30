package sort

// put the configuration methods for sorting here

import (
	"fmt"
	"time"
	"log"
	fastq "squish/fastq"
)

// default values; make sure these match!
const defaultSortMethod string = "alpha"
const defaultSortDescrption string = "Alphabetical sort on sequence"
// for some reason having trouble accessing the const's from main.go...
func GetDefaultSortMethod()(string, string){
	return defaultSortMethod, defaultSortDescrption
}

// holds all the parameters that are used for sorting
type Config struct {
	SortMethod       SortMethod
	InputFilepath    string
	InputFileSize    int64
	OutputFilepath   string
	RecordDelim      byte
	RecordHeaderChar byte
	TimeStart        time.Time
	OrderFilename    string
}

func (config *Config) Run(){
	log.Printf("Using sort method: %v\n", config.SortMethod.CLIArg)
	config.SortMethod.Runner(*config)
}

// keep the CLIarg's associated with sort method and runner methods
type SortMethod struct {
	CLIArg      string
	Description string
	Func        func(*[]fastq.FastqRead)
	Runner func(Config)
}

// parse & return the available sorting methods
func GetSortingMethods() (map[string]SortMethod, string) {
	sortMethodMap := map[string]SortMethod{
		// alpha sort by default!
		defaultSortMethod: SortMethod{
			defaultSortMethod,
			defaultSortDescrption,
			SortReadsSequence,
			RunSortLoader}, // NOTE: MAKE SURE THE DEFAULT LOADER IS CORRECT FOR THE DEFAULT METHOD!
		"gc":              SortMethod{
			"gc", "GC Content Sort", SortReadsGC, RunSortLoader},
		"qual":            SortMethod{
			"qual", "Quality score sort", SortReadsQual, RunSortLoader},
		"alpha-heap":      SortMethod{
			"alpha-heap", "Sequence alpha heap sort", HeapSortSequence, RunSortLoader},
		// "kmer": SortMethod{"kmer", "Kmer sort", _sort.SortKmer},
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
