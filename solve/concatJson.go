package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"

	"github.com/eritema/concatemers/lib"
	"github.com/eritema/concatemers/lib/Alignments"
	"github.com/eritema/concatemers/lib/fasta"
)

type Block struct {
	Id                 int
	ReadLen            int
	BlockLen           int
	BlockDirection     string
	VectorBlockStart   int
	VectorBlockStop    int
	ReadBlockStart     int
	ReadBlockStop      int
	VectorBlockSequnce string
	ReadBlockSequence  string
}

type Read struct {
	Id     string
	Blocks []Block
}

func ClearBlock(block Block) Block {
	block.BlockDirection = ""
	block.BlockLen = 0
	block.ReadLen = 0
	block.VectorBlockStart = 0
	block.VectorBlockStop = 0
	block.ReadBlockStart = 0
	block.VectorBlockStop = 0
	block.VectorBlockSequnce = ""
	block.ReadBlockSequence = ""
	return block
}

func isPointIn(point lib.Point, pointSet []lib.RectSet, strand string) []int {

	var fragments []int
	if strand != "minus" {
		for i := 0; i < len(pointSet); i++ {
			p := pointSet[i]
			flag := lib.InRectPlus(p.PointLists[len(p.PointLists)-1], point, w, d)
			if flag == 1 {
				fragments = append(fragments, i)
			}
		}
	} else {
		for i := 0; i < len(pointSet); i++ {
			p := pointSet[i]
			flag := lib.InRectMinus(p.PointLists[len(p.PointLists)-1], point, w, d)
			if flag == 1 {
				fragments = append(fragments, i)
			}
		}
	}
	return fragments
}

// w is the admitted mismatch number (mm). It should be w=k+mm-1
var wFlag = flag.Int("w", 4, "Number of mismatches alloved (4)")
var dFlag = flag.Int("d", 1, "Number of insertions/deletions alloved (2)")
var kFlag = flag.Int("k", 10, "k-mer dimension (10)")
var thrFlag = flag.Int("thr", 20, "Minimum dimension of a segment (20)")
var directionFlag = flag.Int("D", 3, "Direction search: 1 (3'LAM), 2 (5'LAM), 3 (both)")
var h1Flag = flag.Int("h1", 0, "Prefix dimension (0)")
var h2Flag = flag.Int("h2", 0, "Prefix dimension (0)")
var fastqFlag = flag.Bool("fq", false, "Input file is in fastq format")
var boolFlag = flag.Bool("dot", false, "Return dotPlot cooridantes (false)")
var jFlag = flag.Bool("j", false, "Write the results in json format (false)")
var fileFlag = flag.String("fIn", "test.fq", "Input file")
var vectFlag = flag.String("fVec", "vector.fa", "Vector fasta file")
var helpFlag = flag.Bool("ver", false, "Return the dotrimm version")
var w, d, k, thr, direction, lenThr int
var printDot, fastq bool

func main() {
	flag.Parse()
	head1 := *h1Flag
	head2 := *h2Flag
	k = *kFlag
	d = *dFlag
	w = *wFlag + k
	thr := *thrFlag
	direction := *directionFlag
	printDot = *boolFlag
	fileIn := *fileFlag
	fileVect := *vectFlag
	isJson := *jFlag
	fastq := *fastqFlag
	help := *helpFlag
	var range1, range2 int
	var block Block
	var read Read
	var reads []Read
	//var maxAlignScore int64
	var pointSetPlus, pointSetMinus []lib.RectSet
	var recPlus, recMinus lib.RectSet
	var p lib.Point
	if help {
		fmt.Println("Solve 2015\nVersion 1.0")
		return
	}
	var fastMap = map[string]fasta.ReadType{}

	if fastq {
		fastMap, _ = fasta.GetFastqFile(fileIn, 0, " ")
	} else {
		fastMap, _ = fasta.GetFastaFile(fileIn, 0, "@")
	}
	_, provirus, _ := fasta.GetFastaSequence(fileVect, 0, "@")
	vectorHash := Alignments.GetSeqHash(provirus, k)
	for key2, val2 := range fastMap {
		read.Id = key2
		//fmt.Println(key2)
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
		points := Alignments.KmersBlocksHashed(vectorHash, len(provirus), val2.Sequence, k, range2)
		if printDot {
			fmt.Println(key2)
		}
		//fmt.Println(points)
		for i := 0; i < range1-k+1; i++ {
			for j := 0; j < range2-k+1; j++ {
				if (points[i][j] == 1 || points[i][j] == 3) && (direction == 1 || direction == 3) {
					if printDot {
						fmt.Println(i+1, j+1, 1)
					}
					p.X = i + 1
					p.Y = j + 1
					fragPlus := isPointIn(p, pointSetPlus, "plus")
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
				if points[i][j] >= 2 && direction >= 2 {
					if printDot {
						fmt.Println(i+1, j+1, -1)
					}
					p.X = i + 1
					p.Y = j + 1
					fragMinus := isPointIn(p, pointSetMinus, "minus")
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
		ClearBlock(block)
		id := 0
		for p := range pointSetPlus {
			// The length of the (+) block is 1 bp shorter
			// dim:=pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].X+(k-1)-pointSetPlus[p].PointLists[0].X
			dim := pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].X + k - pointSetPlus[p].PointLists[0].X
			if len(pointSetPlus[p].PointLists) > 0 && dim >= thr {
				block.Id = id
				id++
				block.BlockLen = dim
				if first == 1 {
					first = 0
					//maxAlignScore=int64(0)
				}
				//deletion:=len(provirus)+1 - (pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].X+k)
				//fmt.Println(pointSetPlus[p].PointLists[0].X-1,pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].X+k-1,":",pointSetPlus[p].PointLists[0].Y-1,pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].Y+k-1)
				fragment1 := provirus[pointSetPlus[p].PointLists[0].X-1 : pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].X+k-1]
				fragment2 := ""
				if pointSetPlus[p].PointLists[0].Y-1 <= pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].Y+k-1 {
					fragment2 = val2.Sequence[pointSetPlus[p].PointLists[0].Y-1 : pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].Y+k-1]
				} else {
					fragment2 = val2.Sequence[pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].Y+k-1 : pointSetPlus[p].PointLists[0].Y-1]
				}
				//sequence1,sequence2,alignScore,identity:=Alignments.GlobalAlign(fragment1,fragment2,int64(0),"normal")
				/*if alignScore>maxAlignScore {
					maxAlignScore=alignScore
				}*/
				if !printDot && !isJson {
					fmt.Print(key2, " ")
				}
				block.ReadLen = len(val2.Sequence)
				block.BlockDirection = "+"
				block.VectorBlockStart = pointSetPlus[p].PointLists[0].X
				block.VectorBlockStop = pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].X + k - 1
				block.ReadBlockStart = pointSetPlus[p].PointLists[0].Y
				block.ReadBlockStop = pointSetPlus[p].PointLists[len(pointSetPlus[p].PointLists)-1].Y + k - 1
				block.VectorBlockSequnce = fragment1
				block.ReadBlockSequence = fragment2
				if !isJson {
					fmt.Println( /*,alignScore," ",identity," ",deletion," ",*/
						block.ReadLen,
						block.BlockLen, /*,float32(identity)/float32(dim)*/
						block.BlockDirection,
						block.VectorBlockStart,
						block.VectorBlockStop,
						block.ReadBlockStart,
						block.ReadBlockStop,
						block.VectorBlockSequnce,
						block.ReadBlockSequence) //," ",sequence1," ",sequence2)
				} else {
					read.Blocks = append(read.Blocks, block)
				}
			}
		}
		ClearBlock(block)
		//fmt.Println(maxAlignScore)
		for p := range pointSetMinus {
			dim := pointSetMinus[p].PointLists[len(pointSetMinus[p].PointLists)-1].X + k - pointSetMinus[p].PointLists[0].X
			if len(pointSetMinus[p].PointLists) > 0 && dim >= thr {
				block.Id = id
				id++
				block.BlockLen = dim
				if first == 1 {
					first = 0
				}
				//deletion:=len(provirus)+1 - (pointSetMinus[p].PointLists[0].Y+k)
				fragment1 := provirus[pointSetMinus[p].PointLists[0].X-1 : pointSetMinus[p].PointLists[len(pointSetMinus[p].PointLists)-1].X+k-1]
				fragment2 := lib.RevStringDNA(val2.Sequence[pointSetMinus[p].PointLists[len(pointSetMinus[p].PointLists)-1].Y-1 : pointSetMinus[p].PointLists[0].Y+k-1])
				//sequence1,sequence2,alignScore,identity:=Alignments.GlobalAlign(fragment1,fragment2,int64(0),"normal")
				if !printDot && !isJson {
					fmt.Print(key2, " ")
				}
				block.ReadLen = len(val2.Sequence)
				block.BlockDirection = "-"
				block.VectorBlockStart = pointSetMinus[p].PointLists[0].X
				block.VectorBlockStop = pointSetMinus[p].PointLists[len(pointSetMinus[p].PointLists)-1].X + k
				block.ReadBlockStart = pointSetMinus[p].PointLists[len(pointSetMinus[p].PointLists)-1].Y
				// The end coordinate of the (-) blocks has one bases more
				// block.ReadBlockStop=pointSetMinus[p].PointLists[0].Y+k
				block.ReadBlockStop = pointSetMinus[p].PointLists[0].Y + k - 1
				block.VectorBlockSequnce = fragment1
				block.ReadBlockSequence = fragment2
				if !isJson {
					fmt.Println( /*alignScore," ",identity," ",*/ len(val2.Sequence), block.BlockLen /*,float32(identity)/float32(dim)*/, block.BlockDirection,
						block.VectorBlockStart,
						block.VectorBlockStop,
						block.ReadBlockStart,
						block.ReadBlockStop,
						block.VectorBlockSequnce, block.ReadBlockSequence) //," ",sequence1," ",sequence2)
				} else {
					read.Blocks = append(read.Blocks, block)
				}
			}
		}
		if first == 0 {
			if isJson {
				reads = append(reads, read)
			}
		} else {
			if !isJson {
				fmt.Println(key2, " NO_VECTOR")
			} else {
				block.VectorBlockSequnce = "NO_VECTOR"
				read.Blocks = append(read.Blocks, block)
				reads = append(reads, read)
			}
		}
		pointSetPlus = pointSetPlus[:0]
		pointSetMinus = pointSetMinus[:0]
		read.Blocks = read.Blocks[:0]
	}
	if isJson {
		jsondata, _ := json.Marshal(reads)
		jsonFile, err := os.Create("./data.json")
		if err != nil {
			fmt.Println(err)
		}
		defer jsonFile.Close()
		jsonFile.Write(jsondata)
		jsonFile.Close()
	}
}
