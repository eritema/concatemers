package lib

import (
	//"fmt"
)

func scoreString(str string) int {
	var col = map [string]int {"A":0,"C":0,"G":0,"T":0,}
	for i:=0;i<len(str);i++ {
		col[str[i:i+1]]++
	}
	max:=0
	for _,val:=range col {
		if val>max {
			max=val
		}
	}
	//fmt.Print(str,"->",len(str),"->",max,"->",len(str)-max,":")
	return len(str)-max
}

func Score (motifs []string) int {
	str:=""
	score:=0
	//fmt.Println(motifs)
	for l:=0;l<len(motifs[0]);l++ {
		for i:=0;i<len(motifs);i++ {
			str=str+motifs[i][l:l+1]
		}
		score+=scoreString(str)
		//fmt.Println(score)
		str=""
	}
	return score
}
