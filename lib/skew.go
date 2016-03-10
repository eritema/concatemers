package lib
import ("strings")


func Skew(dna string) ([]int,[]int) {
	dna=strings.ToUpper(dna)	
	var skewMin[]int
	var skew []int
	min:=2
	i:=0
	skew=append(skew,0)
	for l:=0;l<len(dna);l++ {
		if (dna[l:l+1]=="C") {
			i--
			if (i<=min) {
				min=i
			}
		}
		if (dna[l:l+1]=="G") {
			i++
		}
		skew=append(skew,i)
	}
	for i=0;i<len(skew);i++ {
		if (skew[i]==min) {
			skewMin=append(skewMin,i)
		}
	}
	return skew,skewMin
}

