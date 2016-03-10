package lib
import (
	 "strings" 
	 //"fmt"
	 )


/*	Count
	Input: Sequence string, k-mer length
	Output: An array of k length patterns, the frequency and the number
*/ 
func Count(text string, k int) (map[string]int,int,int) {
    c:=0
    count:=0
    var m=map[string]int {}
    text=strings.ToUpper(strings.Replace(text," ","",-1)) // Replace newline chars
    text=strings.ToUpper(strings.Replace(text,"\n","",-1)) // Replace space chars 
    for i:=0;i<len(text)-k+1;i++ {
      m[text[i:i+k]]++
      if m[text[i:i+k]] > c { // the word has a maximal count
        c=m[text[i:i+k]]
        count=1
      } else if m[text[i:i+k]] == c {
      	count++
      }
    }
    return  m,c,count
}

/*	CountMismatches
	Input: Sequence string, k-mer length, number of mismatches
	Output: An array of k length patterns, the frequency and the number
*/
func CountMismatches(text string, k int, d int) (map[string]int,int,int) {
	max:=0
    count:=0
    var n=map[string]int {}
    var m=map[string]int {}
    for i:=0;i<len(text)-k+1;i++ {
    		neighbors:=Neighbors(text[i:i+k],d)
    		for _,keys:=range(neighbors) {
				m[keys]++
			}
    		for keys,_:=range(m) {
					n[keys]++
					if n[keys] > max { // the word has a maximal count
        				max=n[keys]
        				count=1
      				} else if n[keys] == max {
      					count++
      				}
      				delete(m,keys)
			} 	
	}
    return n,max,count
}


/*
	PatternPositions
	Input: Sequence string, pattern string
	Output: an array of positions where the pattern is found
	
*/
func PatternPositions(text string, pattern string) []int {
	var positions []int
	k:=len(pattern)
	//fmt.Println(pattern)
	for i:=0;i<len(text)-k+1;i++ {
		//fmt.Println(text[i:i+k]," ",pattern)
		if string(text[i:i+k])==pattern {
			positions=append(positions,i)
		}
	}
	return positions
}


/*
	PatternPositionsMismatches
	Input: Sequence string, pattern string, number of mismatches
	Output: an array of positions where the pattern is found
	
*/
func PatternPositionsMismatches(text string,pattern string, d int) []int {
	var positions []int
	var m = map[string]int {}
	neighbors:=Neighbors(pattern,d)
	for _,keys:=range(neighbors) {
		m[keys]++
	}
	for keys,_:=range(m) {
		pos:=PatternPositions(text,keys)
		if len(pos)>0 {
			for i:=0;i<len(pos);i++ {
				positions=append(positions,pos[i])
			}
		}
	}
	return positions
}

/*
	CountPatternMismatches
	Input: Sequence string, pattern string, number of mismatches
	Output: Frequency in witch the pattern is found in sequence 
	
*/
func CountPatternMismatches(text string, pattern string, d int) int {
	return len(PatternPositionsMismatches(text,pattern,d))
}

