package main

import (
	"runtime"
	"flag"
	"fmt"
	"os"
	"strings"
	"github.com/eritema/concatemers/lib"
	"github.com/eritema/concatemers/lib/Alignments"
	"github.com/eritema/concatemers/lib/fasta"
)
func doTrimming(head1, head2 int, provirus string, key2 string, val2 fasta.ReadType, pointSetPlus, pointSetMinus []lib.RectSet, recPlus, recMinus lib.RectSet) {
	var p lib.Point

	if head1 == 0 {
		range1 = len(provirus)
	} else {
		range1 = head1
	}
	if head2 == 0 {
		range2 = len(val2.Sequence)
	} else {
		if head2 < len(val2.Sequence) {
			range2 = head2
		} else {
			range2 = len(val2.Sequence)
		}
	}
	points := Alignments.KmersBlocks(provirus, val2.Sequence, k, range1, range2)
	if printDot {
		fmt.Println(key2)
	}
	//fmt.Println(points)
	for i := 0; i < range1-k+1; i++ {
		for j := 0; j < range2-k+1; j++ {
			if (points[i][j] == 1 || points[i][j] == 3) && !direction {
				if printDot {
					fmt.Println(i+1, j+1, 1)
				}
				p.X = i + 1
				p.Y = j + 1
				fragPlus := lib.IsPointIn(p, pointSetPlus, "plus",w,d)
				//fmt.Println(i,j,fragPlus,fragMinus)
				if len(fragPlus) > 0 {
					for fragPlusIndex := 0; fragPlusIndex < len(fragPlus); fragPlusIndex++ {
						pointSetPlus[fragPlus[fragPlusIndex]].PointLists = append(pointSetPlus[fragPlus[fragPlusIndex]].PointLists, p)
					}
				} else {
					recPlus.PointLists = []lib.Point{p}
					pointSetPlus = append(pointSetPlus, recPlus)
				}
			}
			if points[i][j] >= 2 && direction {
				if printDot {
					fmt.Println(i+1, j+1, -1)
				}
				p.X = i + 1
				p.Y = j + 1
				fragMinus := lib.IsPointIn(p, pointSetMinus, "minus",w,d)
				//fmt.Println(p,fragMinus)
				//fmt.Println(i,j,fragPlus,fragMinus)
				if len(fragMinus) > 0 {
					for fragMinusIndex := 0; fragMinusIndex < len(fragMinus); fragMinusIndex++ {
						pointSetMinus[fragMinus[fragMinusIndex]].PointLists = append(pointSetMinus[fragMinus[fragMinusIndex]].PointLists, p)
					}
				} else {
					recMinus.PointLists = []lib.Point{p}
					pointSetMinus = append(pointSetMinus, recMinus)
				}
			}
		}
	}
	if printDot {
		fmt.Println()
	}
	first := 1
	deleTrue := 0
	trimEnd := 0
	if !direction {
		minDel := len(provirus) + 1
		for p := range pointSetPlus {
			dim := pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].X + (k - 1) - pointSetPlus[p].PointLists[0].X
			if len(pointSetPlus[p].PointLists) > 0 && dim >= thr {
				if first == 1 {
					first = 0
				}
				deletion := len(provirus) + 1 - (pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].X + k)
				//fmt.Println(pointSetPlus[p].pointLists[0].X-1,pointSetPlus[p].pointLists[len(pointSetPlus[p].pointLists)-1].X+k-1,":",pointSetPlus[p].pointLists[0].Y-1,pointSetPlus[p].pointLists[len(pointSetPlus[p].pointLists)-1].Y+k-1)
				if deletion < minDel {
					minDel = deletion
					trimEnd = pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].Y + k - 1
				}
			}
		}
		deleTrue = minDel
		pointSetPlus = pointSetPlus[:0]
	} else {
		maxDel := -100000
		for p := range pointSetMinus {
			dim := pointSetMinus[p].PointLists[len(pointSetMinus[p].PointLists)-1].X + k - pointSetMinus[p].PointLists[0].X
			if len(pointSetMinus[p].PointLists) > 0 && dim >= thr {
				if first == 1 {
					first = 0
				}
				deletion := pointSetMinus[p].PointLists[0].X - 1
				if deletion > maxDel {
					maxDel = deletion
					trimEnd = pointSetMinus[p].PointLists[0].Y + k - 1
				}
			}
		}
		deleTrue = maxDel
		pointSetMinus = pointSetMinus[:0]
	}

	if first == 0 {
		if trimEnd < len(val2.Sequence) && len(val2.Sequence[trimEnd:]) >= lenThr {
			if fastq {
				fmt.Printf("%s ReadLen:%v DelBases:%v TrimEnd:%v RemainingBases:%v\n%s\n+\n%s\n", key2, len(val2.Sequence), deleTrue, trimEnd, len(val2.Sequence[trimEnd:]), val2.Sequence[trimEnd:], val2.Quality[trimEnd:])
				//fmt.Println(val2.Sequence)
			} else {
				fmt.Printf("%s ReadLen:%v DelBases:%v TrimEnd:%v RemainingBases:%v\n%s\n", key2, len(val2.Sequence), deleTrue, trimEnd, len(val2.Sequence[trimEnd:]), val2.Sequence[trimEnd:])
			}
		} else {
			totalRemoved++
		}
	} else {
		totalRemoved++
		//fmt.Println(key2," NO_LTR")
	}
	totalReads++
	return
}



// Global variables
var wFlag = flag.Int("w", 1, "Number of mismatches alloved (1)")
var dFlag = flag.Int("d", 1, "Number of insertions/deletions alloved (2)")
var kFlag = flag.Int("k", 10, "k-mer dimension (10)")
var thrFlag = flag.Int("thr", 20, "Minimum dimension of a segment (20)")
var lenThrFlag = flag.Int("l", 20, "Minimum dimension of a valid fasta read (20)")
var itrFlag = flag.Bool("I", false, "The viral terminal repeats are inverted (false)")
var h1Flag = flag.Int("h1", 0, "Vector Prefix dimension (0)")
var h2Flag = flag.Int("h2", 0, "Reads Prefix dimension (0)")
var fastqFlag = flag.Bool("fq", false, "Input file is in fastq format")
var fragmentFlag = flag.String("frag", "", "Terminal DNA sequence to trimm")
var boolFlag = flag.Bool("dot", false, "Return dotPlot cooridantes (false)")
var fileFlag = flag.String("fIn", "test.fq", "Input file")
var vectFlag = flag.String("fVec", "vector.fa", "Vector fasta file")
var helpFlag = flag.Bool("ver", false, "Return the dotrimm version")
var concFlag = flag.Int("p", 1, "Number of concurrent processes (1)")
var w, d, k, thr, lenThr int
var printDot, fastq,direction bool
var totalReads, totalRemoved int
var range1, range2 int



func main() {
	runtime.GOMAXPROCS(runtime.NumCPU())
	// Assigns input control parameters
	flag.Parse()
	head1 := *h1Flag
	head2 := *h2Flag
	k = *kFlag
	d = *dFlag
	w = *wFlag + k
	thr = *thrFlag
	direction = *itrFlag
	printDot = *boolFlag
	fileIn := *fileFlag
	fragment := *fragmentFlag
	fileVect := *vectFlag
	lenThr = *lenThrFlag
	fastq := *fastqFlag
	help := *helpFlag
	conc:= *concFlag


	var pointSetPlus, pointSetMinus []lib.RectSet
	var recPlus, recMinus lib.RectSet
	c := make(chan int,50)
	var provirus string
	if help {
		fmt.Println("dotrim 2015.09.07\nVersion 1.1")
		return
	}
	var fastMap = map[string]fasta.ReadType{}
	if fastq {
		fastMap, _ = fasta.GetFastqFile(fileIn, 0, " ")
	} else {
		fastMap, _ = fasta.GetFastaFile(fileIn, 0, "@")
	}
	if fragment == "" {
		_, provirus, _ = fasta.GetFastaSequence(fileVect, 0, "@")
	} else {
		provirus = strings.ToUpper(fragment)
	}
	if conc<=1 {
		fmt.Fprintf(os.Stderr,"concurrencies %v\n",conc)
		for key2, val2 := range fastMap {
			doTrimming(head1, head2, provirus, key2, val2, pointSetPlus, pointSetMinus,recPlus,recMinus)
		}
	} else {
		fmt.Fprintf(os.Stderr,"concurrencies %v\n",conc)
		for key2, val2 := range fastMap {
			go func() {
				doTrimming(head1, head2, provirus, key2, val2, pointSetPlus, pointSetMinus,recPlus,recMinus)
				c<-1
			}()
			<-c
		}
	}
	fmt.Fprintf(os.Stderr, "Total Reads: %v\tTotal Removed Reads: %v\n", totalReads, totalRemoved)
	//fmt.Println(totalReads,totalRemoved)
}
