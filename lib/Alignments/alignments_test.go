package Alignments

import "testing"

func TestEdit(t *testing.T) {
	seq1 := "ATTACGATAGCAATAGCAGTACTGCAGT"
	seq2 := "ATTACGATAGCAATAGCAGTACTGCAGT"
	seq3 := "ATTACGATAGCAATACCAGTACTGCAGT"
	seq4 := "TTACGATAGCAATACCAGTACTGCAGT"
	if EditDistance(seq1, seq2) != 0 {
		t.Error(`EditDistance("ATTACGATAGCAATAGCAGTACTGCAGT","ATTACGATAGCAATAGCAGTACTGCAGT") not 0`)
	}
	if EditDistance(seq1, seq3) != 1 {
		t.Error(`EditDistance("ATTACGATAGCAATAGCAGTACTGCAGT","ATTACGATAGCAATACCAGTACTGCAGT") not 0`)
	}
	if EditDistance(seq1, seq4) != 2 {
		t.Error(`EditDistance("ATTACGATAGCAATAGCAGTACTGCAGT","TTACGATAGCAATACCAGTACTGCAGT") not 0`)
	}
}
