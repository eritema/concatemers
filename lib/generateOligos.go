package lib

/* This function generate the complete space of k-mer oligos
	recursively
*/

func Space(buque []string,k int) []string {
	alphabet:="ACGT"
	var out []string
	if len(buque)==0 {
		for l:=0;l<4;l++ {
			out=append(out,alphabet[l:l+1])
			
		}
	} else {
		for n:=0;n<len(buque);n++ {
			for l:=0;l<4;l++ {
				out=append(out,buque[n]+alphabet[l:l+1])
				
			}
		}
	}	
	if len(out[0])<k	{	
		return Space(out,k)
	} else {
		return out
	}
}

