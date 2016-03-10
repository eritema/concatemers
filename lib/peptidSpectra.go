package lib
import ("sort";"strconv")

/* Table of the integer mass of 20 standard AA */
var AAMass=map[string]int {
	"G":57,
	"A":71,
	"S":87,
	"P":97,
	"V":99,
	"T":101,
	"C":103,
	"I":113,
	"L":113,
	"N":114,
	"D":115,
	"K":128,
	"Q":128,
	"E":129,
	"M":131,
	"H":137,
	"F":147,
	"R":156,
	"Y":163,
	"W":186,
}
var AAMassRed=map[string]int {
	"57":57,
	"71":71,
	"87":87,
	"97":97,
	"99":99,
	"101":101,
	"103":103,
	"113":113,
	"114":114,
	"115":115,
	"128":128,
	"129":129,
	"131":131,
	"137":137,
	"147":147,
	"156":156,
	"163":163,
	"186":186,
}


type intArray []int
func (s intArray) Len() int { return len(s) }
func (s intArray) Swap(i, j int) { s[i], s[j] = s[j], s[i] }
func (s intArray) Less(i, j int) bool { return s[i] < s[j] }

/*
	SubCyclePeptides return all the substrings from 
	a circular version of the peptide pep
*/
func SubCyclePeptides (pep string) []string {
	var fragment string
	var subPeptides []string
	subPeptides=append(subPeptides,pep)
	for n:=1;n<len(pep);n++ {
		for i:=0;i<len(pep);i++ {
			if(i+n)>len(pep)-1 {
				fragment=pep[i:]+pep[0:(i+n-len(pep))]
			} else {
				fragment=pep[i:i+n]
			}
			subPeptides=append(subPeptides,fragment)
		}
	}
	return subPeptides
}

func CirculPermut(pep string)[]string {
	var fragment string
	var subPeptides []string
	subPeptides=append(subPeptides,pep)
	for i:=0;i<len(pep);i++ {
		fragment=pep[i:]+pep[0:i]
		subPeptides=append(subPeptides,fragment)
	}
	return subPeptides
}

/*
  MolDalton return the molecular weigth of a given peptide
  using the reduced alphabet
*/
func MolDaltonRed(peptide []int) int {
	dalton:=0
	for i:=0;i<len(peptide);i++ {
			dalton+=peptide[i]
	}
	return dalton
}

/*
  MolDalton return the molecular weigth of a given peptide
*/
func MolDalton(peptide string) int {
	dalton:=0
	for i:=0;i<len(peptide);i++ {
			dalton+=AAMass[peptide[i:i+1]]
	}
	return dalton
}

/*
SubSpectrum return true if the integer spectraum spec1 is contained in spec2.
*/
func SubSpectrum (spec1 []int, spec2 []int) bool {
	var found bool
	for i:=0;i<len(spec1);i++ {
		found=false
		for j:=0;j<len(spec2);j++ {
			if spec1[i]==spec2[j] {
					found=true
			}
		}
		if (!found) {
				break
		}	
	}
	return found
}

/* 
Linearspectrum return the theoretical spectrum of the given linear peptide
*/
func Linearspectrum(peptide string) []int {
	var spectra []int
	prefix:=make([]int,len(peptide)+1) 
	prefix[0]=0
	for i:=1;i<len(peptide)+1;i++ {
		prefix[i]=MolDalton(peptide[0:i])
	}
	spectra=append(spectra,0)
	for i:=0;i<len(peptide);i++ {
		for j:=i+1;j<len(peptide)+1;j++ {
			spectra=append(spectra,prefix[j]-prefix[i])
		}
	}	
	sort.Sort(intArray(spectra))
	return spectra
}

/* 
Cyclospectrum return the theoretical spectrum of the given circular peptide
*/
func Cyclospectrum(peptide string) []int {
	var spectra []int
	spectra=append(spectra,0)
	subPeptides:=SubCyclePeptides(peptide)	
	for p:=0;p<len(subPeptides);p++ {
		spectra=append(spectra,MolDalton(subPeptides[p]))
	}
	
	sort.Sort(intArray(spectra))
	return spectra
}

/*
Score return the number of integer masses shared between the a give peptide and integer spectrum
*/


func ScoreIntInt(theoSpectrum []int, spectrum []int) int {
	score:=0
	for i:=0;i<len(theoSpectrum);i++ {
		for j:=0;j<len(spectrum);j++ {
			if theoSpectrum[i]==spectrum[j] {
				spectrum[j]=-1
				score++	
				break
			}
		}
	}
	return score
}

func ScoreStringInt(peptide string, src []int) int {
	var theoSpectrum []int
	spectrum := make([]int, len(src))
	copy(spectrum,src)
	score:=0
	theoSpectrum=Cyclospectrum(peptide)
	for i:=0;i<len(theoSpectrum);i++ {
		for j:=0;j<len(spectrum);j++ {
			if theoSpectrum[i]==spectrum[j] {
				spectrum[j]=-1
				score++	
				break
			}
		}
	}
	return score
}

func ScoreStringIntLinear(peptide string, src []int) int {
	var theoSpectrum []int
	spectrum := make([]int, len(src))
	copy(spectrum,src)
	score:=0
	theoSpectrum=Linearspectrum(peptide)
	for i:=0;i<len(theoSpectrum);i++ {
		for j:=0;j<len(spectrum);j++ {
			if theoSpectrum[i]==spectrum[j] {
				spectrum[j]=-1
				score++	
				break
			}
		}
	}
	return score
}


func AAtoMol (peptide string) string {
	var mol string
	mol=strconv.Itoa(MolDalton(peptide[0:1]))
	for i:=1;i<len(peptide);i++ {
		mol=mol+"-"+strconv.Itoa(MolDalton(peptide[i:i+1]))
	}
	return mol
}