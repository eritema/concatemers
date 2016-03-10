package Alignments

import (
	"strings"
	//"fmt"
 	"github.com/eritema/concatemers/lib"
 	"github.com/eritema/concatemers/lib/strFunctions"
)

var sequenceAGlob=""
var sequenceBGlob=""
var identityGlobal=0

func max(a, b, c int64) (int64,int64) {
    if a >= b {
        if a > c {
        	return a,-1
        }
        return c,0
    } else if b > c {
    	return b,1
    } else {
    	return c,0
   	}
}

func outputLCS(b [][]int64,v string,w string, i int,j int) {
	
	if i==0 {
		for  {
			if j==0 {
				break
			}
			sequenceAGlob=sequenceAGlob+"-"
			sequenceBGlob=sequenceBGlob+w[j-1:j]
			j--
		}
		return
	}
	
	if j==0 {
		for  {
			if i==0 {
				break
			}
			sequenceAGlob=sequenceAGlob+v[i-1:i]
			sequenceBGlob=sequenceBGlob+"-"
			i--
		}
		return
	}
	if b[i][j]==-1 {
		sequenceAGlob=sequenceAGlob+v[i-1:i]
		sequenceBGlob=sequenceBGlob+"-"
		outputLCS(b,v,w,i-1,j)
		
	} else if b[i][j]==1 {	
		sequenceAGlob=sequenceAGlob+"-"
		sequenceBGlob=sequenceBGlob+w[j-1:j]
		outputLCS(b,v,w,i,j-1)
	} else {
		//fmt.Print(v[i-1:i])
		sequenceAGlob=sequenceAGlob+v[i-1:i]
		sequenceBGlob=sequenceBGlob+w[j-1:j]
		if v[i-1:i]==w[j-1:j] {
			identityGlobal++
		}
		outputLCS(b,v,w,i-1,j-1)
	}
}

func GlobalAlign(seq1,seq2 string, sigma int64, blosumString string) (string,string,int64,int)  {
	identityGlobal=0
	var blosum=map[string]int64 {}
	sequenceAGlob=""
	sequenceBGlob=""
	if strings.ToUpper(blosumString)=="BLOSUM62" {
		blosum=lib.GetBlosum("$GOPATH/Data/blosum62.txt")
	} else if strings.ToUpper(blosumString)=="PAM250" {
		blosum=lib.GetBlosum("/media/sf_workspace/BioinfoAlgos/src/github.com/eritema/Data/pam250.txt")
	} else {
		blosum=lib.GetBlosum("/media/sf_workspace/BioinfoAlgos/src/github.com/eritema/Data/simpleScore.txt")
	}
	s := make([][]int64, len(seq1)+1)
	b := make([][]int64, len(seq1)+1)
	for i := range s {
        s[i] = make([]int64, len(seq2)+1)
        b[i] = make([]int64, len(seq2)+1)
    }
	s[0][0]=0
	for i:=1;i<=len(seq1);i++ {
		s[i][0]=s[i-1][0]-sigma
	}
	for j:=1;j<=len(seq2);j++ {
		s[0][j]=s[0][j-1]-sigma
	}
	for i:=1;i<=len(seq1);i++{
		for j:=1;j<=len(seq2);j++ {
			//fmt.Println(s[i-1][j],s[i][j-1],s[i-1][j-1])
			s[i][j],b[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]+blosum[seq1[i-1:i]+seq2[j-1:j]])
			
		}	
	}
	outputLCS(b,seq1,seq2,len(seq1),len(seq2))
	return strConversion.ReverseString(sequenceAGlob),strConversion.ReverseString(sequenceBGlob),s[len(seq1)][len(seq2)],identityGlobal
}

