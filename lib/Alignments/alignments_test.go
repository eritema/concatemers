package Alignments

import "testing"

func TestEdit(t *testing.T) {
	var tests = []struct {
		input1 string
		input2 string
		output int64
	}{
		{"ATTACGATAGCAATAGCAGTACTGCAGT", "ATTACGATAGCAATAGCAGTACTGCAGT", 0},
		{"ATTACGATAGCAATAGCAGTACTGCAGT", "ATTACGATAGCAATACCAGTACTGCAGT", 1},
		{"ATTACGATAGCAATAGCAGTACTGCAGT", "TTACGATAGCAATACCAGTACTGCAGT", 2},
	}
	for _, test := range tests {
		if got := EditDistance(test.input1, test.input2); got != test.output {
			t.Errorf("EditDistance(%q,%q)=%d", test.input1, test.input2, got)
		}
	}
}
