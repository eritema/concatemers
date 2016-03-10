package Alignments

import (
)

var edit int

func minEdit(a, b, c int64, s1 string, s2 string) int64 {
	if s1!=s2 {
		c=c+1
	}	
	a+=1
	b+=1
    if a <= b {
        if a < c {
        	return a
        }
        return c
    } else if b < c {
    	return b
    } else {
    	return c
   	}
}


// Return the edit distance between seq1 and seq 2.
// Mismatches and insertion/deletion have cost of +1
func EditDistance(seq1, seq2 string ) int64 {
	edit=0
	
	s := make([][]int64, len(seq1)+1)
	for i := range s {
        s[i] = make([]int64, len(seq2)+1)
    }
	s[0][0]=0
	for i:=1;i<=len(seq1);i++ {
		s[i][0]=int64(i)
	}
	for j:=1;j<=len(seq2);j++ {
		s[0][j]=int64(j)
	}
	for i:=1;i<=len(seq1);i++{
		for j:=1;j<=len(seq2);j++ {
			s[i][j]=minEdit(s[i-1][j],s[i][j-1],s[i-1][j-1],seq1[i-1:i],seq2[j-1:j])
		}	
	}	
	return s[len(seq1)][len(seq2)]
}