package Alignments

import (
	//"fmt"
	"github.com/eritema/concatemers/lib"
)



// Given a string seq1 and an integer K it returns an hash table of all the k-mers of seq1 

func GetSeqHash(seq1 string, k int) map[string][]int {
	var k1=map[string][]int {}
	for i:=0;i<len(seq1)-k+1;i++ {
		k1[seq1[i:i+k]]=append(k1[seq1[i:i+k]],i)
	}
	return k1
}


// Return a dotplot matrix (int) of the seq1 and seq2 given:
// - the k-mer dimension
// - the length of the fragment in which to serach the kmers (0 all the sequence)

func KmersBlocks(seq1,seq2 string, k int, head1,head2 int) [][]int {
	var k1=map[string][]int {}
	var k2=map[string][]int {}
	var k3=map[string][]int {}
	
	if head1>0 {
		seq1=seq1[:head1]
	}
	if head2>0 {
		seq2=seq2[:head2]
	}
	rev:=lib.RevStringDNA(seq2)
	points := make([][]int, len(seq1))
	for i := range points {
        points[i] = make([]int, len(seq2))
    }
	for i:=0;i<len(seq1)-k+1;i++ {
		k1[seq1[i:i+k]]=append(k1[seq1[i:i+k]],i)
	}
	for i:=0;i<len(seq2)-k+1;i++ {
		k2[seq2[i:i+k]]=append(k2[seq2[i:i+k]],i)
		k3[rev[i:i+k]]=append(k3[rev[i:i+k]],len(rev)-i-k)
	}
	
	for key,val1:=range(k1) {
		if val2,ok:=k2[key]; ok {
			for i:=0;i<len(val1);i++ {
				for j:=0;j<len(val2);j++{
					points[val1[i]][val2[j]]+=1
				}
			}
		}
		if val3,ok:=k3[key]; ok {
			for i:=0;i<len(val1);i++ {
				for j:=0;j<len(val3);j++{
					points[val1[i]][val3[j]]+=2
					//fmt.Println("in")
				}
			}
		}
	}
	return points
}

// Return a dotplot matrix (int) of k1 (hast k-mers table) and seq2 (string) given:
// - the k-mer dimension
// - the length of the fragment in which to serach the kmers (0 all the sequence)

func KmersBlocksHashed(k1 map[string][]int,vectLen int,seq2 string, k int, head2 int) [][]int {
	var k2=map[string][]int {}
	var k3=map[string][]int {}
	
	if head2>0 {
		seq2=seq2[:head2]
	}
	rev:=lib.RevStringDNA(seq2)
	points := make([][]int, vectLen)
	for i := range points {
        points[i] = make([]int, len(seq2))
    }
	for i:=0;i<len(seq2)-k+1;i++ {
		k2[seq2[i:i+k]]=append(k2[seq2[i:i+k]],i)
		k3[rev[i:i+k]]=append(k3[rev[i:i+k]],len(rev)-i-k)
	}
	
	for key,val1:=range(k1) {
		if val2,ok:=k2[key]; ok {
			for i:=0;i<len(val1);i++ {
				for j:=0;j<len(val2);j++{
					points[val1[i]][val2[j]]+=1
				}
			}
		}
		if val3,ok:=k3[key]; ok {
			for i:=0;i<len(val1);i++ {
				for j:=0;j<len(val3);j++{
					points[val1[i]][val3[j]]+=2
					//fmt.Println("in")
				}
			}
		}
	}
	return points
}
