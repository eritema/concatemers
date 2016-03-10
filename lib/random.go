package lib

import (
   "math/rand"
   //"fmt"
)


func Random (t int32) int {
	
	return int(rand.Int31n(t))
}

func RandomVect(prob []float64) int {
	var i int
	var cumul []float64
	cumul=append(cumul,0.0)
	j:=1
	var sum float64
	for i=0;i<len(prob);i++ {
		cumul=append(cumul,cumul[j-1]+prob[i])
		sum+=prob[i]
		j++
	}
	for i=1;i<len(cumul);i++ {
		cumul[i]=cumul[i]/sum
	}
	pos:=rand.Float64()
	//fmt.Print(cumul," ",pos)
	for i=0;i<len(cumul);i++ {
		if pos>cumul[i] && pos<cumul[i+1] {
			return i
		}
	}
	return i
}