package lib

import (

)

var profIndex = map [string]int{
		"A":0,
		"C":1,
		"G":2,
		"T":3,
	}

func GetMostProbMotif(profile [4][]float64,sequence string,k int) string {
	var probable=""
	var max=0.0
	prob:=1.0
	for i:=0;i<len(sequence)-k+1;i++ {
		for l:=0;l<k;l++ {
			prob*=profile[profIndex[sequence[i+l:i+l+1]]][l]
		}
		if (prob>max) {
			max=prob
			probable=sequence[i:i+k]
		}
		prob=1.0
	}
	if probable=="" {
		return sequence[0:k]
	}
	return probable
}

func GetProbMotifProb(profile [4][]float64,sequence string, k int) []float64 {
	var probable []float64
	prob:=1.0
	for i:=0;i<len(sequence)-k+1;i++ {
		for l:=0;l<k;l++ {
			prob*=profile[profIndex[sequence[i+l:i+l+1]]][l]
		}
		probable=append(probable,prob)
		prob=1.0
	}
	return probable
}