/*
 * hikeeba-personafi.go
 *
 * Based on 'simplegrep.go' example.
 *
 * Use hyperscan to fuzzy regex match flanking sequences
 * for subsequent processing.
 *
 * Brett Whitty <bwhitty@pgdx.com>
 *
 */
package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"regexp"
	_ "runtime"
	"strconv"
	"strings"

	//hyperscan
	"github.com/flier/gohs/hyperscan"

	//Heng Li's FASTQ file reader
	//"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/drio/drio.go/bio/fasta"

	// logging
	log "github.com/sirupsen/logrus"

	// progress bar DEBUG
	"github.com/cheggaaa/pb/v3"
	_ "time"
	//^^^^^^^^^^^^^^^^

	_ "github.com/davecgh/go-spew/spew" //debugging
)

var (
	// to be set during build
	Binary    string
	Version   string
	BuildDate string
	DebugFlag string
	Debug     bool

	// flags
	flagNoColor    = flag.Bool("C", false, "Disable colorized output.")
	flagRecompile  = flag.Bool("c", false, "Force pattern database recompile.")
	flagPrintId    = flag.Bool("i", false, "Print ID of sequence record.")
	flagFASTQOut   = flag.Bool("q", false, "Print FASTQ output.")
	flagFASTQMSeq  = flag.Bool("m", false, "Include matched sequence in FASTQ / ID output formats.")
	flagLTrim      = flag.Bool("L", false, "Trim sequence left/upstream of match.")
	flagRTrim      = flag.Bool("R", false, "Trim sequence right/downstream of match.")
	flagMTrim      = flag.Bool("M", false, "Trim matched sequence.")
	flagRevComp    = flag.Bool("r", false, "Reverse-complement output.")
	flagByteOffset = flag.Bool("b", false, "Display offset in bytes of a matched pattern.")
	flagProgress   = flag.Bool("progress", true, "Show progress bar.")
	flagSilent     = flag.Bool("s", false, "Silent mode.")

	fileWriter []io.Writer
)

var theme = func(s string) string { return s }

func init() {
	// setting DebugFlag = false will cause parameters
	// with values that don't pass validation to be deleted
	Debug, _ = strconv.ParseBool(DebugFlag)

	// init logger
	log.SetOutput(os.Stderr)

	if Debug {
		log.SetLevel(log.DebugLevel)
	}

	flag.Parse()
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

// Record contains the data from a fasta fastq record
type FASTQRecord struct {
	Name, Seq, Qual string
}

//type FASTQRecord interface {
//}

/**
 * This is the function that will be called for each match that occurs. @a ctx
 * is to allow you to have some application-specific state that you will get
 * access to for each match. In our simple example we're just going to use it
 * to pass in the pattern that was being searched for so we can print it out.
 */
func eventHandler(id uint, from, to uint64, flags uint, context interface{}) error {
	// TODO: can maybe be optimized
	inputData := []byte(strings.TrimSpace(context.(FASTQRecord).Seq) + "\n")
	inputQual := []byte(strings.TrimSpace(context.(FASTQRecord).Qual) + "\n")

	start := bytes.LastIndexByte(inputData[:from], '\n')
	end := int(to) + bytes.IndexByte(inputData[to:], '\n')

	if start == -1 {
		start = 0
	} else {
		start += 1
	}

	if end == -1 {
		end = len(inputData)
	}

	if *flagByteOffset {
		fmt.Printf("%d: ", start)
	}

	seqLeftString := string(inputData[:from])
	seqMatchString := string(inputData[from:to])
	seqRightString := string(inputData[to:end])
	qualLeftString := string(inputQual[:from])
	qualMatchString := string(inputQual[from:to])
	qualRightString := string(inputQual[to:end])

	if *flagLTrim {
		seqLeftString = ""
		qualLeftString = ""
	}
	if *flagRTrim {
		seqRightString = ""
		qualRightString = ""
	}
	if *flagMTrim {
		seqMatchString = ""
		qualMatchString = ""
	}

	// support modes later
	outSeq := seqLeftString + seqMatchString + seqRightString
	outQual := qualLeftString + qualMatchString + qualRightString
	outName := context.(FASTQRecord).Name

	reId := regexp.MustCompile(`^@(\S+)`)
	outId := reId.FindStringSubmatch(outName)[1]

	if *flagRevComp {
		outSeq = reverseComplementDNA(outSeq)
		outQual = reverse(outQual)
	}

	if *flagFASTQOut {
		matchSeq := ""
		if *flagFASTQMSeq {
			matchSeq = " " + seqMatchString
		}
		//fmt.Printf("%s\n%s\n%s\n%s\n", outName, outSeq, "+"+outId+" "+fmt.Sprint(id)+":"+seqMatchString, outQual)
		fmt.Printf("%s\n%s\n%s\n%s\n", outName, outSeq, "+"+outId+" "+fmt.Sprint(id)+":"+fmt.Sprint(from)+"-"+fmt.Sprint(to)+matchSeq, outQual)
	} else if *flagPrintId {
		matchSeq := ""
		if *flagFASTQMSeq {
			matchSeq = " " + seqMatchString
		}
		fmt.Printf("%s\n", outId+" "+fmt.Sprint(id)+":"+fmt.Sprint(from)+"-"+fmt.Sprint(to)+matchSeq)
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
		fmt.Fprint(os.Stderr, "HIKEEBA! ["+fmt.Sprintf("%s %s(%s) DEBUG=%t", Binary, Version, BuildDate, Debug)+"] // Brett Whitty <bwhitty@pgdx.com>\n")
	}

	if flag.NArg() != 2 {
		fmt.Fprintf(os.Stderr, "Usage: %s ["+green("flags")+"] <"+cyan("pattern file")+"> <"+cyan("input file")+">\n", highlight(os.Args[0]))
		os.Exit(-1)
	}

	if !*flagNoColor {
		stat, _ := os.Stdout.Stat()

		if stat != nil && stat.Mode()&os.ModeType != 0 {
			theme = highlight
		}
	}

	//pattern := hyperscan.NewPattern(flag.Arg(0), hyperscan.DotAll|hyperscan.SomLeftMost)
	//pattern := hyperscan.NewPattern(flag.Arg(0), hyperscan.SomLeftMost|hyperscan.Caseless)
	patternFile := flag.Arg(0)

	//    spew.Dump(hyperscan.ExpressionInfo(pattern.Expression))
	//   return

	// input file
	inputFN := flag.Arg(1)

	// input is STDIN?
	isSTDIN := strings.EqualFold("/dev/stdin", inputFN) || strings.EqualFold("stdin", inputFN) || strings.EqualFold("-", inputFN)

	// check if input is gzipped
	isGzip, err := regexp.MatchString(`\.gz$`, inputFN)
	checkErr(err)

	// file size in bytes
	var fileSizeB int64 = 0

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
	fmt.Fprintf(os.Stderr, highlight("INFO: ")+"Pattern file: %s\n", patternFile)
	//dbStreaming, dbBlock := databasesFromFile(patternFile)
	database := blockDatabaseFromFile(patternFile)

	defer database.Close()

	var inFile *os.File
	// open file for reading
	if isSTDIN {
		// input is STDIN
		inFile = os.Stdin
	} else {
		fileInfo, err := os.Stat(inputFN)
		checkErr(err)
		fileSizeB = fileInfo.Size()

		// input is normal file
		inFile, err = os.Open(inputFN)
		checkErr(err)
	}
	defer inFile.Close()

	var bar *pb.ProgressBar
	if fileSizeB != 0 {
		bar = pb.Full.Start64(fileSizeB)
	}

	// set up a FASTQ reader, supporting gz'd stream
	var fqr fasta.FqReader
	if isGzip {
		barReader := bar.NewProxyReader(inFile)
		gzipReader, err := gzip.NewReader(barReader)
		checkErr(err)
		fqr.Reader = bufio.NewReader(gzipReader)
	} else {
		// get buffered reader
		barReader := bar.NewProxyReader(inFile)
		fqr.Reader = bufio.NewReader(barReader)
	}
	defer bar.Finish()

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

	seqCount := 0
	for r, done := fqr.Iter(); !done; r, done = fqr.Iter() {
		seqCount++

		// will pass this into the eventHandler
		record := FASTQRecord{Name: r.Name, Seq: r.Seq, Qual: r.Qual}

		//go func() {

		// create byte buffer from FASTQ record sequence string
		// => strings.TrimSpace() may be overkill here
		// eventHandler is expecting input is a line terminated with "\n"
		inputData := []byte(strings.TrimSpace(record.Seq) + "\n")

		//		fmt.Fprintf(os.Stderr, "%d: %s\n", lineCount, inputData)
		//fmt.Fprintf(os.Stderr, "Scanning %d bytes with Hyperscan\n", len(inputData))

		if err := database.Scan(inputData, scratch, eventHandler, record); err != nil {
			fmt.Fprint(os.Stderr, "ERROR: Unable to scan input buffer. Exiting.\n")
			os.Exit(-1)
		}

		//}()
	}

	/* Scanning is complete, any matches have been handled, so now we just
	 * clean up and exit.
	 */

	return
}

func scanFastqRecord(database hyperscan.BlockDatabase, scratch *hyperscan.Scratch, record FASTQRecord) {

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
		lineno += 1

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
		fmt.Fprintf(os.Stderr, highlight("INFO: ")+"Serialized pattern DB exists, use '-c' flag to recompile from text.\n")
		fmt.Fprintf(os.Stderr, highlight("INFO: ")+"Reading from pattern DB file: %s\n", dbFilename)
		return readDbFile(dbFilename)
	}

	fmt.Fprintf(os.Stderr, highlight("INFO: ")+"Compiling patterns ... ")

	// do the actual file reading and pattern parsing
	patterns := parseFile(filename)

	// compile a block database from the parsed patterns
	bdb, err := hyperscan.NewBlockDatabase(patterns...)
	checkErr(err, fmt.Sprintf("Could not compile patterns, %s", err))

	fmt.Fprintf(os.Stderr, "DONE!\n")

	// serialize the compiled database to save time on next run
	fmt.Fprintf(os.Stderr, highlight("INFO: ")+"Serializing pattern DB ... ")
	writeDbFile(dbFilename, bdb)
	fmt.Fprintf(os.Stderr, "DONE!\n")

	return bdb
}

// returns string to use as serialized pattern database file name
func getDbFilename(filename string) string {
	return filename + ".hsdb"
}

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

	// open a filehandle for writing the DB file
	outFile, err := os.Create(dbFilename)
	checkErr(err, fmt.Sprintf("Couldn't open DB file '%s' for writing! %s", dbFilename, err))
	defer outFile.Close()

	// create gzip stream writer
	gzWriter, err := gzip.NewWriterLevel(outFile, gzip.BestCompression)
	checkErr(err, fmt.Sprintf("Couldn't create gzip writer! %s", err))
	defer gzWriter.Close()

	// write gzip compressed DB bytes to file
	gzWriter.Write(dbData)
}

func getGzWriter(filename string) io.Writer {
	// open a filehandle for writing the file
	outFile, err := os.Create(filename)
	checkErr(err, fmt.Sprintf("Couldn't open file '%s' for writing! %s", outFile, err))

	// create gzip stream writer
	gzWriter, err := gzip.NewWriterLevel(outFile, gzip.BestCompression)
	checkErr(err, fmt.Sprintf("Couldn't create gzip writer! %s", err))

	return gzWriter
}
