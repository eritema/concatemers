package main

import (
	"encoding/json"
	"fmt"
	//"os"
	"flag"
	"io/ioutil"
	"log"
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

// w is the admitted mismatch number (mm). It should be w=k+mm-1
var fileFlag = flag.String("fIn", "test.fq", "Input file")
var expandFlag = flag.Bool("e", false, "Show the vector regions in the fragments")
var thrLenFlag = flag.Int("t", 20, "Minimum lenght of a synthenic block")
var helpFlag = flag.Bool("ver", false, "Return the dotrimm version")

func main() {
	regionsMap := map[int]string{
		0: "trs",
		1: "a",
		2: "C",
		3: "3A",
		4: "c",
		5: "b",
		6: "RBEp",
		7: "B",
		8: "A",
	}

	regions := [][]int{
		{0, 18, 32},
		{1, 39, 62},
		{2, 63, 71},
		{3, 72, 74},
		{4, 75, 83},
		{5, 85, 87},
		{6, 0, 0},
		{7, 92, 94},
		{8, 95, 140},
	}

	flag.Parse()
	fileIn := *fileFlag
	help := *helpFlag
	thrLen := *thrLenFlag
	expand := *expandFlag
	var data []Read
	//var maxAlignScore int64
	if help {
		fmt.Println("reconstruct 2015\nVersion 1.0")
		return
	}
	file, err := ioutil.ReadFile(fileIn)
	if err != nil {
		log.Fatal(err)
		return
	}
	err = json.Unmarshal(file, &data)
	if err != nil {
		log.Fatal(err)
		return
	}
	for i := 0; i < len(data); i++ {
		fmt.Println(data[i].Id)
		helpBlocks := make([]Block, len(data[i].Blocks))
		for j := 0; j < len(data[i].Blocks); j++ {
			helpBlocks[data[j].Blocks.BlockReadStart] = data[j].Blocks
		}
		for j := 0; j < len(helpBlocks); j++ {
			primo := true

			//fmt.Println(data[i].Blocks[j].Id,data[i].Blocks[j].VectorBlockStart,data[i].Blocks[j].VectorBlockStop)
			//fmt.Println(data[i].Blocks[j].Id)
			if data[i].Blocks[j].BlockDirection == "+" && data[i].Blocks[j].BlockLen >= thrLen {
				fmt.Print(data[i].Blocks[j].BlockLen, "\t", data[i].Blocks[j].BlockDirection, "\t", data[i].Blocks[j].ReadBlockStart, "\t", data[i].Blocks[j].ReadBlockStop, "\t")
				for l := 0; l < len(regions); l++ {
					if data[i].Blocks[j].VectorBlockStart > regions[l][2] {
						continue
					}
					if primo && data[i].Blocks[j].VectorBlockStart < regions[l][2] && data[i].Blocks[j].VectorBlockStop > regions[l][2] {

						fmt.Print(data[i].Blocks[j].VectorBlockStart)
						primo = false
					}
					if data[i].Blocks[j].VectorBlockStop >= regions[l][2] {
						if expand {
							fmt.Print(regionsMap[l], ":")
						}
						continue
					} else {
						if expand {
							fmt.Print(regionsMap[l], data[i].Blocks[j].VectorBlockStop)
							fmt.Println()
						} else {
							fmt.Print(":", data[i].Blocks[j].VectorBlockStop)
							fmt.Println()
						}
						break
						//						fmt.Println(data[i].Blocks[j].Id,"->",data[i].Blocks[j].BlockLen)
					}
				}
			} else if data[i].Blocks[j].BlockDirection == "-" && data[i].Blocks[j].BlockLen >= thrLen {
				fmt.Print(data[i].Blocks[j].BlockLen, "\t", data[i].Blocks[j].BlockDirection, "\t", data[i].Blocks[j].ReadBlockStart, "\t", data[i].Blocks[j].ReadBlockStop, "\t")
				for l := 0; l < len(regions); l++ {
					if data[i].Blocks[j].VectorBlockStart > regions[l][2] {
						continue
					}
					if primo && data[i].Blocks[j].VectorBlockStart < regions[l][2] && data[i].Blocks[j].VectorBlockStop > regions[l][2] {

						fmt.Print(data[i].Blocks[j].VectorBlockStart)
						primo = false
					}
					if data[i].Blocks[j].VectorBlockStop >= regions[l][2] {
						if expand {
							fmt.Print(regionsMap[l], ":")
						}
						continue
					} else {
						if expand {
							fmt.Print(regionsMap[l], data[i].Blocks[j].VectorBlockStop)
							fmt.Println()
						} else {
							fmt.Print(":", data[i].Blocks[j].VectorBlockStop)
							fmt.Println()
						}
						break
						//						fmt.Println(data[i].Blocks[j].Id,"->",data[i].Blocks[j].BlockLen)
					}
				}
			}
		}
	}
}
