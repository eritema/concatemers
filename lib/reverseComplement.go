package lib

// complemetarity table for DNA
var revDNA=map[string]string {
	"a":"T",
	"A":"T",
	"c":"G",
	"C":"G",
	"g":"C",
	"G":"C",
	"t":"A",
	"T":"A",
	"N":"N",
	"n":"N"}
// complemetarity table for RNA
var revRNA=map[string]string {
	"a":"U",
	"A":"U",
	"c":"G",
	"C":"G",
	"g":"C",
	"G":"C",
	"u":"A",
	"U":"A",
	"N":"N",
	"n":"N"}

// Return the revers complement of a file
func RevFile(buffer []byte) string {
	ReversedString:=""
	s:=string(buffer[:])
	for count:=len(s)-1;count>=0;count-- {
		ReversedString+=revDNA[string(s[count:count+1])]
	}
	return ReversedString	
}

// Return the reverse complement of a DNA sequence
func RevStringDNA(buffer string) string {
	ReversedString:=""
	for count:=len(buffer)-1;count>=0;count-- {
		ReversedString+=revDNA[string(buffer[count:count+1])]
	}
	return ReversedString	
}

// Return the reverse complement of a RNA sequence
func RevStringRNA(buffer string) string {
	ReversedString:=""
	for count:=len(buffer)-1;count>=0;count-- {
		ReversedString+=revRNA[string(buffer[count:count+1])]
	}
	return ReversedString	
}