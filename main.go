/*
 *   Copyright (c) 2020
 *   All rights reserved.
 */

package main

/*
 * hikeeba.go
 *
 * HIKEEBA! GoBCLy
 * => FASTQ "floating barcodes" demultiplexing functionality
 *
 * ( TODO: factor this out into a 'hikeeba gobcly fastq' sub-command using cobra )
 *
 * 2019-20, Brett Whitty <bwhitty@pgdx.com>
 *
 */

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"crypto/md5"
	"encoding/hex"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	//hyperscan
	//Heng Li's FASTQ file reader => not using
	// logging
	// progress bar
	// DEBUG: not used
	_ "runtime"
	_ "time"
	//^^^^^^^^^^^^^^^^
	//

	"github.com/flier/gohs/hyperscan"
	// TODO: re-evaluate FASTQ readers vs. line reading
	//"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/cheggaaa/pb/v3"
	_ "github.com/davecgh/go-spew/spew" //debugging
	"github.com/drio/drio.go/bio/fasta"
	_ "github.com/gobuffalo/packr"
	log "github.com/sirupsen/logrus"
)

var (
	// Binary = Compile-time name of binary
	Binary    string
	// Cmd = Compile-time name of command; TODO: TEMP for debug
	Cmd       string = "GoBCLy"
	// Version = Compile-time version string
	Version   string            
	// BuildDate = Compile-time date string
	BuildDate string
	// DebugFlag = Compile-time debug string
	DebugFlag string
	              // Debug flag
	Debug     bool

	// flags

	// input flags
	flagR1File       = flag.String("r1", "", "Path to R1 file.")
	flagR2File       = flag.String("r2", "", "Path to R2 file.")
	flagPatternsFile = flag.String("p", "", "Path to Hyperscan-complatible PCRE patterns table file.")

	// database flags
	flagRecompile = flag.Bool("c", false, "Force pattern database recompile.")

	// global output options
	// => match / non-match output options
	flagLTrim   = flag.Bool("L", false, "Trim sequence left/upstream of match.")
	flagRTrim   = flag.Bool("R", false, "Trim sequence right/downstream of match.")
	flagMTrim   = flag.Bool("M", false, "Trim matched sequence.")
	flagRevComp = flag.Bool("r", false, "Reverse-complement output.")
	// => FASTQ output options
	flagFASTQOut  = flag.Bool("q", false, "Print FASTQ output.")
	flagFASTQMSeq = flag.Bool("m", true, "Include matched sequence in FASTQ / ID output formats.")
	// generic match output
	flagPrintID    = flag.Bool("i", false, "Print ID of sequence record.")
//	flagByteOffset = flag.Bool("b", false, "Display offset in bytes of a matched pattern.")

	// logging / debug options
	flagNoColor = flag.Bool("C", false, "Disable colorized output.")
	flagDebug   = flag.Bool("d", false, "Debug mode.")
	flagSilent  = flag.Bool("s", false, "Silent mode.")
	
	// fileWriters map of GzipWriters
	fileWriters GzipWriters
)
// InputSet = Set of required files for validation purposes
type InputSet struct {
	R1Filepath   string
	R2Filepath   string
	PatternsFile string
	OK           bool
} // TODO: this isn't necessary, remove later

	// GzipWriters for storing pointers to io.Writers
	type GzipWriters map[uint]*gzip.Writer
	// => index type 'uint' here matches 'id' type of hyperscan match

	var theme = func(s string) string { return s }

func init() {
	fileWriters = make(GzipWriters)

//	box := packr.NewBox("./.packr")

	flag.Parse()

	// setting DebugFlag = false will cause parameters
	// with values that don't pass validation to be deleted
	Debug, _ = strconv.ParseBool(DebugFlag)

	// in case compile-time replacement doesn't happen
	if Binary == "" {
		Binary = os.Args[0]
		Version = "dev"
		BuildDate = "1976-10-19"
		Debug = true
		DebugFlag = "true"
	}

	if *flagDebug {
		Debug = true
	}
	if !*flagSilent {
		Debug = false
	}

	// init logger
	log.SetFormatter(&log.TextFormatter{
		DisableColors: false,
		FullTimestamp: true,
	})
	log.SetLevel(log.InfoLevel)
	log.SetOutput(os.Stderr)

	// TODO: flag for log redirect to file

	if Debug {
		log.SetLevel(log.DebugLevel)
	} else if *flagSilent {
		log.SetLevel(log.ErrorLevel)
	}
}

func cyan(s string) string {
	return "\033[36m" + s + "\033[0m"
	//    return s
}

func green(s string) string {
	return "\033[32m" + s + "\033[0m"
	//    return s
}
func highlight(s string) string {
	return "\033[35m" + s + "\033[0m"
	//    return s
}

// FASTQRecord contains the data from a fasta fastq record
type FASTQRecord struct {
	InputFileBasename, Name, Seq, Qual string
}

// DemuxRecord probably isn't necessary
type DemuxRecord struct {
	R1 FASTQRecord
	R2 FASTQRecord
} // TODO: remove this

// DemuxReaders also probably aren't needed
type DemuxReaders struct {
	R1 *fasta.FqReader
	R2 *fasta.FqReader
} //TODO: remove this

// DemuxScratch hyperscan scratch space
type DemuxScratch struct {
	R1 *hyperscan.Scratch
	R2 *hyperscan.Scratch
} //TODO: remove this

var demuxSet [4]FASTQRecord

//type FASTQRecord interface {
//}

/**
 * This is the function that will be called for each match that occurs. @a ctx
 * is to allow you to have some application-specific state that you will get
 * access to for each match. In our simple example we're just going to use it
 * to pass in the pattern that was being searched for so we can print it out.
 */
func eventHandler(id uint, from, to uint64, flags uint, context interface{}) error {

	//context.(FASTQRecord).Seq
	fastq := FASTQRecord{context.(FASTQRecord).InputFileBasename, context.(FASTQRecord).Name, context.(FASTQRecord).Seq, context.(FASTQRecord).Qual}
	//fastq := FASTQRecord{InputFileBasename: "R1", Name: fqR1.Name, Seq: fqR1.Seq, Qual: fqR1.Qual}
	
	// TODO: can maybe be optimized
	//inputBasename := []byte(strings.TrimSpace(context.(FASTQRecord).FileBasename) + "\n")
	inputData := []byte(strings.TrimSpace(fastq.Seq) + "\n")
	inputQual := []byte(strings.TrimSpace(fastq.Qual) + "\n")

	// determine match start and end positions
	matchStartPos := bytes.LastIndexByte(inputData[:from], '\n')
	// fix start bounds
	if matchStartPos == -1 {
		matchStartPos = 0
	} else {
		matchStartPos++
	}
	matchEndPos := int(to) + bytes.IndexByte(inputData[to:], '\n')
	// fix end bounds
	if matchEndPos == -1 {
		matchEndPos = len(inputData)
	}

	/*
	TODO: remove or find some use for this
	
	if *flagByteOffset {
		fmt.Printf("%d: ", matchStartPos)
	}
	*/

	var seqLeftString, qualLeftString, seqMatchString, qualMatchString, seqRightString, qualRightString string

	// optionally trim sequence left / upstream of match
	if *flagLTrim {
		seqLeftString = ""
		qualLeftString = ""
	} else {
		seqLeftString = string(inputData[:from])
		qualLeftString = string(inputQual[:from])
	}
	// optionally trim sequence that was matched
	// TODO: masking options
	if *flagMTrim {
		seqMatchString = ""
		qualMatchString = ""
	} else {
		seqMatchString = string(inputData[from:to])
		qualMatchString = string(inputQual[from:to])
	}
	// optionally trim sequence right / downstream of match
	if *flagRTrim {
		seqRightString = ""
		qualRightString = ""
	} else {
		seqRightString = string(inputData[to:matchEndPos])
		qualRightString = string(inputQual[to:matchEndPos])
	}
	
	// prepare output strings for writing
	// TODO: mode specific?
	outSeq := seqLeftString + seqMatchString + seqRightString
	outQual := qualLeftString + qualMatchString + qualRightString
	outName := fastq.Name

	// reverse complement the output if user requested
	if *flagRevComp {
		outSeq = reverseComplementDNA(outSeq)
		outQual = reverse(outQual)
	}

	// use a regex to parse out FASTQ read ID
	// TODO: Fix this hack; this is FASTQ specific
	reID := regexp.MustCompile(`^@(\S+)`)
	outID := reID.FindStringSubmatch(outName)[1]

	if *flagFASTQOut {
		if fileWriters[id] == nil {
			outputGzFastqFile := fastq.InputFileBasename + "." + fmt.Sprintf("%d", id) + ".hs_dmux.fastq.gz"
			fileWriters[id] = getGzWriter(outputGzFastqFile)
		}
		gzWriter := fileWriters[id]

		matchSeq := ""
		if *flagFASTQMSeq {
			matchSeq = " " + seqMatchString
		}
		//fmt.Printf("%s\n%s\n%s\n%s\n", outName, outSeq, "+"+outID+" "+fmt.Sprint(id)+":"+seqMatchString, outQual)
		gzWriter.Write([]byte(fmt.Sprintf("%s\n%s\n%s\n%s\n", outName, outSeq, "+"+outID+" "+fmt.Sprint(id)+":"+fmt.Sprint(from)+"-"+fmt.Sprint(to)+matchSeq, outQual)))
	} else if *flagPrintID {
		matchSeq := ""
		if *flagFASTQMSeq {
			matchSeq = " " + seqMatchString
		}
		fmt.Printf("%s %s\n", fastq.InputFileBasename, outID+" "+fmt.Sprint(id)+":"+fmt.Sprint(from)+"-"+fmt.Sprint(to)+matchSeq)
	} else {

		// DEBUG fmt.Printf("start=%d, end=%d, from=%d, to=%d\n", start, end, from, to)
		// DEBUG fmt.Printf("%s%s%s\n", inputData[start:from], theme(string(inputData[from:to])), inputData[to:end])
		//		fmt.Printf("%s%s\n", theme(string(inputData[from:to])), inputData[to:end])
		fmt.Printf("%s%s%s\n", seqLeftString, theme(seqMatchString), seqRightString)
	}

	return nil
}

func main() {
	/*
		// TEST CODE FOR BUILDING DATABASE FROM MULTIPLE EXPRESSIONS
		b := hyperscan.DatabaseBuilder{}
		db, err := b.AddExpressions("101:/abc/Q", "102:/def/Q", "/(101&102)/Co").Build()
		info, err := db.Info()
		mode, err := info.Mode()
		spew.Dump(b, db, info, mode)
		db.Close()
		// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	*/

	//panic(fmt.Sprintf("PANIC!\n"))

	if !*flagSilent {
		fmt.Fprint(os.Stderr, highlight("HIKEEBA!")+" "+cyan(Cmd)+" "+"["+fmt.Sprintf("%s %s(%s) DEBUG=%t", Binary, Version, BuildDate, Debug)+"] // Brett Whitty <bwhitty@pgdx.com>\n")
	}

	if *flagR1File == "" || *flagR2File == "" {
		fmt.Fprintf(os.Stderr, "Usage: %s ["+green("flags")+"] <"+cyan("pattern file")+"> <"+cyan("input file")+">\n", highlight(Binary))
		flag.PrintDefaults()
		os.Exit(-1)
	}

	if !*flagNoColor {
		// enable color if supported, unless disabled by flag
		stat, _ := os.Stdout.Stat()
		if stat != nil && stat.Mode()&os.ModeType != 0 {
			theme = highlight
		}
	} else {
		log.SetFormatter(&log.TextFormatter{
			DisableColors: true,
		})
	}

	//pattern := hyperscan.NewPattern(flag.Arg(0), hyperscan.DotAll|hyperscan.SomLeftMost)
	//pattern := hyperscan.NewPattern(flag.Arg(0), hyperscan.SomLeftMost|hyperscan.Caseless)
	patternFile := *flagPatternsFile

	//    spew.Dump(hyperscan.ExpressionInfo(pattern.Expression))
	//   return

	// input file
	// TODO: legacy
	//inputFN := *flagR1File

	inputSet := InputSet{*flagR1File, *flagR2File, *flagPatternsFile, false}

	// TODO:	log.

	// TODO: new
	//	inputR1File := *flagR1File
	//	inputI1File := *flagI1File
	//	inputR2File := *flagR2File
	//	inputI2File := *flagI2File

	//	var barR1, barI1, barR2, barI2 int64
	//	var bar *pb.ProgressBar

	//	var dFQR DemuxReaders
	readerR1, bar := getFQReader(inputSet.R1Filepath, true)
	readerR2, _ := getFQReader(inputSet.R2Filepath, false)

	// TODO: fix this hack
	r1FileBasename, r1OK := getGzFastqBasename(inputSet.R1Filepath)
	if !r1OK {
		log.Fatal("R1 file doesn't have '.fastq.gz' suffix as expected!")
	}
	r2FileBasename, r2OK := getGzFastqBasename(inputSet.R2Filepath)
	if !r2OK {
		log.Fatal("R2 file doesn't have '.fastq.gz' suffix as expected!")
	}

	//	totalBytes := bytesR1 + bytesI1 + bytesR2 + bytesI2
	//	bar = pb.Full.Start64(totalBytes)

	// TODO: factor this out into 'getFileType'

	// input is STDIN?
	//	isSTDIN := strings.EqualFold("/dev/stdin", inputFN) || strings.EqualFold("stdin", inputFN) || strings.EqualFold("-", inputFN)

	// check if input is gzipped
	//	isGzip, err := regexp.MatchString(`\.gz$`, inputFN)
	//	checkErr(err)

	// TODO: move this into some sort of file integrity check?
	// file size in bytes
	//	var inputFileSizeBytes int64 = 0

	// TODO: ADD EDIT / HAMMING DISTANCE SUPPORT!!!

	//pattern.HammingDistance := 1

	/* First, we attempt to compile the pattern provided on the command line.
	 * We assume 'DOTALL' semantics, meaning that the '.' meta-character will
	 * match newline characters. The compiler will analyse the given pattern and
	 * either return a compiled Hyperscan database, or an error message
	 * explaining why the pattern didn't compile.
	 */
	//database, err := hyperscan.NewBlockDatabase(pattern) // <= compile pattern
	//if err != nil {
	//	fmt.Fprintf(os.Stderr, "ERROR: Unable to compile pattern \"%s\": %s\n", pattern.String(), err.Error())
	//	os.Exit(-1)
	//}

	// Read our pattern set in and build Hyperscan databases from it.
	log.Info(fmt.Sprintf("Pattern file: %s\n", patternFile))
	//dbStreaming, dbBlock := databasesFromFile(patternFile)
	database := blockDatabaseFromFile(patternFile)
	defer database.Close()

	/* if false {
		var inFile *os.File
		// open file for reading
		if isSTDIN {
			// input is STDIN
			inFile = os.Stdin
		} else {
			inputFileSizeBytes = getFileSizeInBytes(inputFN)

			// input is normal file
			inFile, err = os.Open(inputFN)
			checkErr(err)
		}
		defer inFile.Close()

		var fileReader io.Reader

		if *flagSilent {
			fileReader = inFile
		} else {
			bar = pb.Full.Start64(inputFileSizeBytes)
			//inFile = bar.NewProxyReader(inFile).(io.Reader)
			fileReader = bar.NewProxyReader(inFile)
			defer bar.Finish()
		}

		// set up a FASTQ reader, supporting gz'd stream
		var fqr fasta.FqReader
		if isGzip {
			//		barReader := bar.NewProxyReader(inFile)
			gzipReader, err := gzip.NewReader(fileReader)
			checkErr(err)
			fqr.Reader = bufio.NewReader(gzipReader)
		} else {
			// get buffered reader
			barReader := bar.NewProxyReader(fileReader)
			fqr.Reader = bufio.NewReader(barReader)
		}

	}
	// buffer to store the sequence data for scanning
	//var inputData []byte

	/* Finally, we issue a call to hs_scan, which will search the input buffer
	 * for the pattern represented in the bytecode. Note that in order to do
	 * this, scratch space needs to be allocated with the hs_alloc_scratch
	 * function. In typical usage, you would reuse this scratch space for many
	 * calls to hs_scan, but as we're only doing one, we'll be allocating it
	 * and deallocating it as soon as our matching is done.
	 *
	 * When matches occur, the specified callback function (eventHandler in
	 * this file) will be called. Note that although it is reminiscent of
	 * asynchronous APIs, Hyperscan operates synchronously: all matches will be
	 * found, and all callbacks issued, *before* hs_scan returns.
	 *
	 * In this example, we provide the input pattern as the context pointer so
	 * that the callback is able to print out the pattern that matched on each
	 * match event.
	*/
	scratch, err := hyperscan.NewScratch(database)
	checkErr(err, fmt.Sprintf("Unable to allocate scratch space. Exiting."))
	defer scratch.Free()

	var demuxScratch DemuxScratch

	// scratch clones for demux threads
	demuxScratch.R1, err = scratch.Clone()
	checkErr(err)
	demuxScratch.R2, err = scratch.Clone()
	checkErr(err)

	doneCount := 0
	seqCountR1 := 0
	seqCountR2 := 0
	for {
		/*
			R1
		*/
		fqR1, doneR1 := readerR1.Iter()
		if doneR1 {
			doneCount++
		}
		rR1 := FASTQRecord{InputFileBasename: r1FileBasename, Name: fqR1.Name, Seq: fqR1.Seq, Qual: fqR1.Qual}
		log.Debug(rR1.Name)
		scanFastqRecord(database, demuxScratch.R1, rR1)
		seqCountR1++
		/*
			R2
		*/
		fqR2, doneR2 := readerR2.Iter()
		if doneR2 {
			doneCount++
		}
		rR2 := FASTQRecord{InputFileBasename: r2FileBasename, Name: fqR2.Name, Seq: fqR2.Seq, Qual: fqR2.Qual}
		log.Debug(rR2.Name)
		scanFastqRecord(database, demuxScratch.R2, rR2)
		seqCountR2++
		if doneCount >= 1 {
			if doneCount == 2 {
				bar.Finish()

				log.Debug(fmt.Sprintf("%d %d", seqCountR1, seqCountR2))

				// close any open gzip filewriters
				for _, fw := range fileWriters {
					fw.Close()
				}

				break
			} else {
				log.Error(fmt.Sprintf("%d %d", seqCountR1, seqCountR2))
				log.Fatal("Encountered input file record mismatch!!!")
				os.Exit(-1)
			}
		}

		// will pass this into the eventHandler
		//		recordR1 := FASTQRecord{Name: rI1.Name, Seq: r.Seq, Qual: r.Qual}

		//go func() {

		// create byte buffer from FASTQ record sequence string
		// => strings.TrimSpace() may be overkill here
		// eventHandler is expecting input is a line terminated with "\n"
		//inputData := []byte(strings.TrimSpace(record.Seq) + "\n")

		// prepare demux recordset for output
		//recordSet := DemuxRecord{record, record, record, record}
		//		recordSet := DemuxRecord{rR1, rI1, rR2, rI2}

		//		for _, record := range []FASTQRecord{recordSet.R1, recordSet.I1, recordSet.R2, recordSet.I2} {

		//			scanFastqRecord(database, scratch, record)
		//		fmt.Fprintf(os.Stderr, "%d: %s\n", lineCount, inputData)
		//fmt.Fprintf(os.Stderr, "Scanning %d bytes with Hyperscan\n", len(inputData))

		///		if err := database.Scan(inputData, scratch, eventHandler, record); err != nil {
		//			fmt.Fprint(os.Stderr, "ERROR: Unable to scan input buffer. Exiting.\n")
		//			os.Exit(-1)
		//		}

		//		}

		//}()
	}

	/* Scanning is complete, any matches have been handled, so now we just
	 * clean up and exit.
	 */

	return
}

func scanFastqRecord(database hyperscan.BlockDatabase, scratch *hyperscan.Scratch, record FASTQRecord) {

	//var demuxRecordKeys = []string("R1","I1","R2","I2")

	//	for i, k := range demuxRecordKeys {
	//		record := demuxRecord[i]

	// => strings.TrimSpace() may be overkill here
	// eventHandler is expecting input is a line terminated with "\n"
	inputData := []byte(strings.TrimSpace(record.Seq) + "\n")

	//		fmt.Fprintf(os.Stderr, "%d: %s\n", lineCount, inputData)
	//fmt.Fprintf(os.Stderr, "Scanning %d bytes with Hyperscan\n", len(inputData))

	if err := database.Scan(inputData, scratch, eventHandler, record); err != nil {
		log.Fatal("ERROR: Unable to scan input buffer. Exiting.")
		os.Exit(-1)
	}

}

//}

func checkErr(err error, optMsg ...string) {
	msg := ""
	if len(optMsg) > 0 {
		msg = optMsg[0]
	}
	if err != nil {
		if msg != "" {
			log.Fatal(msg)
		}
		os.Exit(-1)
	}
}

func reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

func reverseComplementDNA(s string) string {
	fwd := "atgcn.ATGCN."
	rvs := "tacgn.TACGN."

	comp := map[rune]rune{
		[]rune(fwd)[0]:  []rune(rvs)[0],
		[]rune(fwd)[1]:  []rune(rvs)[1],
		[]rune(fwd)[2]:  []rune(rvs)[2],
		[]rune(fwd)[3]:  []rune(rvs)[3],
		[]rune(fwd)[4]:  []rune(rvs)[4],
		[]rune(fwd)[5]:  []rune(rvs)[5],
		[]rune(fwd)[6]:  []rune(rvs)[6],
		[]rune(fwd)[7]:  []rune(rvs)[7],
		[]rune(fwd)[8]:  []rune(rvs)[8],
		[]rune(fwd)[9]:  []rune(rvs)[9],
		[]rune(fwd)[10]: []rune(rvs)[10],
		[]rune(fwd)[11]: []rune(rvs)[11],
	}

	runes := []rune(s)
	// check that we have a mapping for all characters we encounter
	for _, c := range runes {
		_, cExists := comp[c]
		if !cExists {
			log.Fatal(fmt.Sprintf("Unexpected character in:\n%s\n", s))
		}
	}

	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = comp[runes[j]], comp[runes[i]]
	}

	return string(runes)
}

func parseFile(filename string) (patterns []*hyperscan.Pattern) {
	// open pattern file for reading
	data, err := ioutil.ReadFile(filename)
	checkErr(err, fmt.Sprintf("Can't read pattern file '%s'", filename))

	reader := bufio.NewReader(bytes.NewBuffer(data))
	eof := false
	lineno := 0

	for !eof {
		line, err := reader.ReadString('\n')

		switch err {
		case nil:
			// pass
		case io.EOF:
			eof = true
		default:
			fmt.Fprintf(os.Stderr, highlight("ERROR: ")+"Can't read pattern file %s, %s\n", filename, err)
			os.Exit(-1)
		}

		line = strings.TrimSpace(line)
		lineno++

		// if line is empty, or a comment, we can skip it
		if len(line) == 0 || line[0] == '#' {
			continue
		}

		// otherwise, it should be ID:PCRE, e.g.
		//  10001:/foobar/is
		strs := strings.SplitN(line, ":", 2)
		// TODO: Need to add support for parsing extended attribute flags...
		//  10001:/foobar/is{key1=value1,key2=value2,...}

		// parse the pattern ID
		id, err := strconv.ParseInt(strs[0], 10, 64)
		checkErr(err, fmt.Sprintf("Could not parse id at line %d, %s", lineno, err))

		// parse the pattern
		pattern, err := hyperscan.ParsePattern(strs[1])
		checkErr(err, fmt.Sprintf("Could not parse pattern at line %d, %s", lineno, err))

		// set Id on the pattern
		pattern.Id = int(id)

		// add the pattern
		patterns = append(patterns, pattern)
	}

	return
}

/**
 * This function will read in the file with the specified name, with an
 * expression per line, ignoring lines starting with '#' and build a Hyperscan
 * database for it.
 */
func blockDatabaseFromFile(filename string) hyperscan.BlockDatabase {
	dbFilename := getDbFilename(filename)

	// remove existing compiled database if recompile flag set
	if *flagRecompile {
		if fileExists(dbFilename) {
			err := os.Remove(dbFilename)
			checkErr(err, fmt.Sprintf("Couldn't remove DB file '%s'! %s", dbFilename, err))
		}
	}

	// short circuit compiling if we already have a serialized db
	if fileExists(dbFilename) {
		log.Debug("Serialized pattern DB exists, use '-c' flag to recompile from text.")
		log.Info(fmt.Sprintf("Reading from pattern DB file: %s", dbFilename))
		return readDbFile(dbFilename)
	}

	log.Info("Compiling patterns ... ")

	// do the actual file reading and pattern parsing
	patterns := parseFile(filename)

	// compile a block database from the parsed patterns
	bdb, err := hyperscan.NewBlockDatabase(patterns...)
	checkErr(err, fmt.Sprintf("Could not compile patterns, %s", err))

	log.Info(" ... DONE!")

	// serialize the compiled database to save time on next run
	log.Info("Serializing pattern DB ... ")
	writeDbFile(dbFilename, bdb)
	log.Info(" ... DONE!")

	return bdb
}

// returns string to use as serialized pattern database file name
func getDbFilename(filename string) string {
	fileMD5, err := getFileMD5(filename)
	checkErr(err, fmt.Sprintf("Failed to generate MD5 digest of %s", filename))
	return filename + "." + fileMD5 + ".hsdb"
}

// check if a file exists by name
func fileExists(filename string) bool {
	info, err := os.Stat(filename)

	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

// reads a pattern block database that we've previously serialized to a gzipped file
func readDbFile(dbFilename string) hyperscan.BlockDatabase {
	// open DB file for reading; should be gzipped
	inFile, err := os.Open(dbFilename)
	checkErr(err, fmt.Sprintf("Couldn't open DB file '%s' for reading! %s", dbFilename, err))
	defer inFile.Close()

	// create gzip stream reader
	gzipReader, err := gzip.NewReader(inFile)
	checkErr(err, fmt.Sprintf("Couldn't create gzip reader on DB file! %s", err))
	defer gzipReader.Close()

	dbData, err := ioutil.ReadAll(gzipReader)
	checkErr(err, fmt.Sprintf("Could not read database file, %s", err))

	bdb, err := hyperscan.UnmarshalBlockDatabase(dbData)
	checkErr(err, fmt.Sprintf("Failed to unmarshall DB file! %s", err))

	return bdb
}

// write a compiled block database to a gzipped file to avoid recompiling on next run
func writeDbFile(dbFilename string, database hyperscan.BlockDatabase) {
	// serialize the block database to bytes
	dbData, err := database.Marshal()
	checkErr(err, fmt.Sprintf("Could not serialize pattern DB, %s", err))

	// get a gzip stream writer
	gzWriter := getGzWriter(dbFilename)
	defer gzWriter.Close()

	// write gzip compressed DB bytes to file
	gzWriter.Write(dbData)
}

// returns a pointer to a gzip.Writer
func getGzWriter(filename string) *gzip.Writer {
	// open a filehandle for writing the file
	outFile, err := os.Create(filename)
	checkErr(err, fmt.Sprintf("Couldn't open file '%s' for writing! %s", filename, err))

	// create gzip stream writer
	gzWriter, err := gzip.NewWriterLevel(outFile, gzip.BestCompression)
	checkErr(err, fmt.Sprintf("Couldn't create gzip writer! %s", err))

	return gzWriter
}

// returns an IO.Reader for a file
func getFileReader(filename string) io.Reader {
	var fileReader io.Reader

	// input is normal file
	fileReader, err := os.Open(filename)
	checkErr(err)

	return fileReader
}

// returns a gzip IO.Reader

// returns the size of a file in bytes
func getFileSizeInBytes(filename string) int64 {
	fileInfo, err := os.Stat(filename)
	checkErr(err)
	fileSizeInBytes := fileInfo.Size()
	return fileSizeInBytes
}

// returns pointer to a FASTQ reader;
// if wantBar = true returns pointer to a pb.ProgressBar for read IO, otherwise nil
func getFQReader(filename string, wantBar bool) (*fasta.FqReader, *pb.ProgressBar) {
	// if wantBar = true, this will be progress bar for read IO
	// ... otherwise = nil
	var bar *pb.ProgressBar

	// input is STDIN?
	isSTDIN := strings.EqualFold("/dev/stdin", filename) || strings.EqualFold("stdin", filename) || strings.EqualFold("-", filename)

	// check if input is gzipped
	isGzip, err := regexp.MatchString(`\.gz$`, filename)
	checkErr(err)

	// TODO: move this into some sort of file integrity check?
	// file size in bytes
	var inputFileSizeBytes int64 = 0

	var inFile *os.File
	// open file for reading
	if isSTDIN {
		// input is STDIN
		inFile = os.Stdin
	} else {
		inputFileSizeBytes = getFileSizeInBytes(filename)

		// input is normal file
		inFile, err = os.Open(filename)
		checkErr(err)
	}

	var fileReader io.Reader

	if *flagSilent {
		fileReader = inFile
	} else {
		if wantBar {
			bar = pb.Full.Start64(inputFileSizeBytes)
			fileReader = bar.NewProxyReader(inFile)
		} else {
			fileReader = inFile
		}
	}

	// set up a FASTQ reader, supporting gz'd stream
	var fqr fasta.FqReader
	if isGzip {
		gzipReader, err := gzip.NewReader(fileReader)
		checkErr(err)
		fqr.Reader = bufio.NewReader(gzipReader)
	} else {
		// get buffered reader
		barReader := bar.NewProxyReader(fileReader)
		fqr.Reader = bufio.NewReader(barReader)
	}

	return &fqr, bar
}

// cut and paste md5 checksum code
func getFileMD5(filePath string) (string, error) {
	var fileDigest string
	file, err := os.Open(filePath)
	if err != nil {
		return fileDigest, err
	}
	defer file.Close()
	hash := md5.New()
	if _, err := io.Copy(hash, file); err != nil {
		return fileDigest, err
	}
	hashInBytes := hash.Sum(nil)[:16]
	fileDigest = hex.EncodeToString(hashInBytes)
	return fileDigest, nil
}

// returns basename of a "*.fastq.gz" file,
// or error if the file doesn't match
func getGzFastqBasename (filePath string) (string, bool) {
	// TODO: do this in a better way so it doesn't fail with caps (eg: regex)
	base := strings.TrimSuffix(filePath, ".fastq.gz")
	if base == filePath {
		return "", false 
	}

	return filepath.Base(base), true
}