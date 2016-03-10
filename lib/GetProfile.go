package lib

import (
//"fmt"
)

var bases = map [string]int {
	"A":0,"C":0,"G":0,"T":0,
}

func GetProfile(motifs []string) [4][]float64 {
	var profile [4][]float64
	
	
	for i:=0;i<len(motifs[0]);i++ {
		for n:=0;n<len(motifs);n++ {
			bases[motifs[n][i:i+1]]++	
		}
		profile[0]=append(profile[0],float64(bases["A"])/float64(len(motifs)))
		bases["A"]=0
		profile[1]=append(profile[1],float64(bases["C"])/float64(len(motifs)))
		bases["C"]=0
		profile[2]=append(profile[2],float64(bases["G"])/float64(len(motifs)))
		bases["G"]=0
		profile[3]=append(profile[3],float64(bases["T"])/float64(len(motifs)))
		bases["T"]=0
	}
	
	return profile
}

func GetProfilePseudo(motifs []string) [4][]float64 {
	var profile [4][]float64
	bases["A"]=1
	bases["C"]=1
	bases["G"]=1
	bases["T"]=1
	for i:=0;i<len(motifs[0]);i++ {
		for n:=0;n<len(motifs);n++ {
			bases[motifs[n][i:i+1]]++
				
		}
		profile[0]=append(profile[0],float64(bases["A"])/float64(len(motifs)+4.0))
		//fmt.Print((bases["A"]),":",len(motifs)+4.0,":",float64(bases["A"])/float64((len(motifs)+4.0))," ")
		bases["A"]=1
		profile[1]=append(profile[1],float64(bases["C"])/float64(len(motifs)+4.0))
		bases["C"]=1
		profile[2]=append(profile[2],float64(bases["G"])/float64(len(motifs)+4.0))
		bases["G"]=1
		profile[3]=append(profile[3],float64(bases["T"])/float64(len(motifs)+4.0))
		bases["T"]=1
	}
	
	return profile
}

func GetProfilePseudoIndex(motifs []string, k int) [4][]float64 {
	var profile [4][]float64
	bases["A"]=1
	bases["C"]=1
	bases["G"]=1
	bases["T"]=1
	for i:=0;i<len(motifs[0]);i++ {
		for n:=0;n<len(motifs);n++ {
			if n==k {continue}
			bases[motifs[n][i:i+1]]++
		}
		profile[0]=append(profile[0],float64(bases["A"])/float64(len(motifs)+4.0))
		//fmt.Print((bases["A"]),":",len(motifs)+4.0,":",float64(bases["A"])/float64((len(motifs)+4.0))," ")
		bases["A"]=1
		profile[1]=append(profile[1],float64(bases["C"])/float64(len(motifs)+4.0))
		bases["C"]=1
		profile[2]=append(profile[2],float64(bases["G"])/float64(len(motifs)+4.0))
		bases["G"]=1
		profile[3]=append(profile[3],float64(bases["T"])/float64(len(motifs)+4.0))
		bases["T"]=1
	}
	
	return profile
}
