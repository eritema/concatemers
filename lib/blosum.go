package lib

import (
	"github.com/eritema/concatemers/lib/strFunctions"
)

// 
func GetBlosum(fileName string) map[string]int64 {
	var blosum=map[string]int64 {}
	
	rows:=GetDataFromFile(fileName)
	aa:=strConversion.GetVectorOfStrings(rows[0], " ")
	for i:=0;i<len(aa);i++ {
		scores:=strConversion.GetVectorOfInt(rows[i+1], " ")
		for j:=0;j<len(aa);j++ {
			blosum[aa[i]+aa[j]]=scores[j]
		}
	}
	return blosum
}	
