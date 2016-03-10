package lib

import ("bytes")


/*
	complementar bases to a certain base
*/
var letters = map[string]string{
	"A":"CGT",
	"C":"AGT",
	"G":"ACT",
	"T":"ACG",
}

var alphabet = "ACGT"

/* ImmediateNeighbors return all the neighbors of the sequence pattern
   that have only one mismatches
   Input: pattern
   Output: array of neighbours of pattern
*/
func ImmediateNeighbors(pattern string) []string {
	var buffer bytes.Buffer
	var Neighbors []string
	for i:=0;i<len(pattern);i++ {
		symbols:=(letters[pattern[i:i+1]])
		for j:=0;j<len(symbols);j++ {
				buffer.WriteString(pattern[0:i])
				buffer.WriteString(symbols[j:j+1])
				buffer.WriteString(pattern[i+1:])
				Neighbors=append(Neighbors,buffer.String())
				buffer.Reset()
		}
	}
	return Neighbors
}

/*	Neighbors return all the neighbours of pattern with max d mismatches
	Input: sequence pattern and number of mismatches d
	Output: array of neigbours
*/


func Neighbors(pattern string, d int) []string {
	var buffer bytes.Buffer
	var neighbors []string
	neighbors=append(neighbors,pattern)
	if (d==0) {return neighbors} 
	if (len(pattern)==1) {
		neighbors=append(neighbors,"A")
		neighbors=append(neighbors,"C")
		neighbors=append(neighbors,"G")
		neighbors=append(neighbors,"T")
		return neighbors 
	}
	//symbols:=(letters[pattern[0:1]])
	suffix:=string(pattern[1:])
	suffixNeighbor:=Neighbors(suffix,d)
	for _,text:=range(suffixNeighbor) {
		hamming,_:=Hamming(text,suffix)
		if(hamming<d) {
			//for j:=0;j<len(symbols);j++ {
			for j:=0;j<len(alphabet);j++ {
				//buffer.WriteString(symbols[j:j+1])
				buffer.WriteString(alphabet[j:j+1])
				buffer.WriteString(text)
				neighbors=append(neighbors,buffer.String())
				buffer.Reset()
			}
		} else {
			buffer.WriteString(pattern[0:1])
			buffer.WriteString(text)
			neighbors=append(neighbors,buffer.String())
			buffer.Reset()
		}
	}
	return neighbors
}