package lib

import (
	"bytes"
	//"fmt"
)
var CodonsRNA = map[string]string{
	"AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T",
	"AGA":"R", "AGC":"S", "AGG":"R", "AGU":"S", "AUA":"I", "AUC":"I", "AUG":"M", "AUU":"I",
	"CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P",
	"CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R", "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L",
	"GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A",
	"GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G", "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V", 
	"UAA":"-", "UAC":"Y", "UAG":"-", "UAU":"Y", "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S", 
	"UGA":"-", "UGC":"C", "UGG":"W", "UGU":"C", "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F",
}
var CodonsDNA = map[string]string{
	"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
	"AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I",
	"CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
	"CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
	"GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
	"GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", 
	"TAA":"-", "TAC":"Y", "TAG":"-", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", 
	"TGA":"-", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F",
}
var CodonsNum = map[string]int {
	"A":4,"L":6,
	"R":6,"K":2,
	"N":2,"M":1,
	"D":2,"F":2,
	"C":2,"P":4,
	"Q":2,"S":6,
	"E":2,"T":4,
	"G":4,"W":1,
	"H":2,"Y":2,
	"I":3,"V":4,
}

/*
TranslateRNA return the 6 frame translation of the provided RNA sequence
*/
func TranslateRNA(rna string) []string {
	var frame bytes.Buffer
	var peptide []string
	for i:=0;i<3;i++ {
		for count:=i;count<len(rna)-3-i;count+=3 {
			frame.WriteString(CodonsRNA[rna[count:count+3]])
		}
		peptide=append(peptide,frame.String())
		frame.Reset()
	}
	rnaRev:=RevStringRNA(rna)
	for i:=0;i<3;i++ {
		for count:=i;count<len(rnaRev)-3-i;count+=3 {
			frame.WriteString(CodonsRNA[rnaRev[count:count+3]])
		}
		peptide=append(peptide,frame.String())
		frame.Reset()
	}
	return peptide
}
/*
	Given a standard 20 AA peptide pep the function return an array of 
	all the possible originating DNA sequences that codify for that pep
	Be aware that the number of fragments grows between 2^N and 4^N, where N
	is the lenght of the peptide
*/
func RetroTranslate(pep string) []string {
	var sequences []string
	var fragment string
	fragment=""
	numSeq:=1
	for i:=0;i<len(pep);i++{
			numSeq=numSeq*CodonsNum[pep[i:i+1]]
	}
	for j:=0;j<numSeq;j++ {
				sequences=append(sequences,fragment)
	}
	precCodNum:=1;
	for i:=0;i<len(pep);i++{
		codNum:=CodonsNum[pep[i:i+1]]
		for n:=0;n<codNum;n++ {
			for j:=0;j<precCodNum;j++ {
				sequences[precCodNum*n+j]=sequences[j]
			}		
		}
		count:=0
		for codon,aa:=range(CodonsDNA){
			if(pep[i:i+1]==aa) {
				for j:=0;j<precCodNum;j++ {
					sequences[precCodNum*count+j]=sequences[precCodNum*count+j]+codon
				}
				count++	
			}
		}
		precCodNum*=codNum	
	}
	return sequences
}

