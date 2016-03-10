package lib

import ("fmt")


func MotifEnumeration(dna map[string]int , k int, d int) map[string]int {
	var motif = map[string]int {}
	var found int
	//per ogni sequenza key in dna
	for key,_:=range(dna) {
		fmt.Println(key)
		//per ogni pattern in key
		for i:=0;i<len(key)-k+1;i++ {
			
			//calcola le stringe a distanza d
			patterns:=Neighbors(key[i:i+k],d)
			
			//per ogni neighbour
			for j:=0;j<len(patterns);j++ {
				
				//vedi se nelle altre key c'e' il pattern
				for keyAltre,_:=range(dna) {
					found=0
					if(key!=keyAltre) {
						//calcola tutti i sottopattern
						for l:=0;l<len(keyAltre)-k+1;l++ {
							ham,_:=Hamming(patterns[j],keyAltre[l:l+k])
							if ham <=d {
								found=1
								break
							}
						}
						if found==0 {break}
					}
				}
				if found==1 {
					motif[patterns[j]]++
				}
			}
		}
	}
	return motif
} 