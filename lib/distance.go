package lib
//import ("fmt")

func Distance(pattern string , dna map[string]int) int {
	dist:=0
	for sequence,_:=range(dna) {
		hamming:=100000
		k:=len(pattern)
		for n:=0;n<len(sequence)-k+1;n++ {
			ham,_:=Hamming(pattern,sequence[n:n+k])
			if ham < hamming {
				hamming=ham
			}
		}
		dist+=hamming
	}	
	return dist
}

func MedianString(dna map[string]int, k int) map[string]int {
	var patterns = []string {}
	var motif = map [string]int {}
	// costruisci lo spazio dei pattern
	minDistance:=10000
	/*for key,_:=range(dna) {
		//riempi lo spazio dei pattern
		for n:=0;n<len(key)-k+1;n++ {
			patterns[key[n:n+k]]++
		}
	}*/
	patterns=Space(patterns,k)
	for i:=0;i<len(patterns);i++ {
			distance:=Distance(patterns[i],dna)
			//fmt.Println(patterns[i]," ",distance)
			if distance<minDistance {
				motif=make(map[string]int)
				motif[patterns[i]]=distance
				minDistance=distance
			} else if distance==minDistance { 
				motif[patterns[i]]=distance
			}	
	}
	return motif
}
