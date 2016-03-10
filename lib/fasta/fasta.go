package fasta

import (
   "github.com/eritema/concatemers/lib"
   "github.com/eritema/concatemers/lib/strFunctions"
   "strings"
   //"fmt"
)

type ReadType struct {
	Sequence string
	Quality string
}

// Read a fasta sequence from file and return a string with header and one with sequence
// if head is 0 the whole header is used othervise a substring delimeited by sep of the header int is used
func GetFastaSequence(fastaFile string,head int, sep string) (string,string,int) {
	var name string
	sequence:=""
	rows:=lib.GetDataFromFile(fastaFile)
	id:=strConversion.GetVectorOfStrings(rows[0], sep)
	if head>0 {
		id[0]=id[0][:head]
	}	
	if id[0][0:1]!=">" {
		name=id[0]
	} else {
		name=id[0][1:]
	}
	for i:=1;i<len(rows); i++{
		sequence=sequence+strings.ToUpper(rows[i])	
	}
	return name,sequence,0
}	


// Read a fasta file and return an dictionary of it
//if head is 0 the whole header is used othervise a substring delimeited by sep of the header int is used
func GetFastaFile(fastaFile string,head int, sep string) (map[string]ReadType,int) {
	var readsDict=map [string]ReadType {}
	var read ReadType
	firstFastaChar:=""
	rows:=lib.GetDataFromFile(fastaFile)
	if head==0 {
		for i:=0;i<len(rows); {
			ids:=strConversion.GetVectorOfStrings(rows[i], sep)
			id:=ids[0]
			if id[0:1]!=">" {
				firstFastaChar=">"
			}
			i++
			read.Sequence=strings.ToUpper(rows[i])
			read.Quality=""
			readsDict[firstFastaChar+id]=read
			i++
		}
		return readsDict,0
	} else {
		for i:=0;i<len(rows); {
			id:=rows[i]
			if id[0:1]!=">" {
				firstFastaChar=">"
			}
			i++
			read.Sequence=strings.ToUpper(rows[i])
			read.Quality=""
			readsDict[firstFastaChar+id[:head]]=read
			i++
		}
		return readsDict,0
	}
}


func GetFastqFile(fastaFile string,head int, sep string) (map[string]ReadType,int) {
	var readsDict=map [string]ReadType {}
	var read ReadType
	firstFastaChar:=""
	rows:=lib.GetDataFromFile(fastaFile)
	if head==0 {
		for i:=0;i<len(rows); {
			ids:=strConversion.GetVectorOfStrings(rows[i], sep)
			//fmt.Println(ids[0])
			id:=ids[0]
			i++
			read.Sequence=strings.ToUpper(rows[i])
			i=i+2
			read.Quality=rows[i]
			readsDict[id]=read
			i++
		}
		return readsDict,0
	} else {
		for i:=0;i<len(rows); {
			id:=rows[i]
			if id[0:1]!="@" {
				firstFastaChar="@"
			}
			i++
			read.Sequence=strings.ToUpper(rows[i])
			i=i+2
			read.Quality=rows[i]
			readsDict[firstFastaChar+id[:head]]=read
			i++
		}
		return readsDict,0
	}
}
