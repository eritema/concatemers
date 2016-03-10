package lib

import (
	"os";"log"
	"bufio";"strings"
	"fmt";"io";"errors"
)



func GetDataFromFile(name string) []string {
	var Dna []string
	file,err := os.Open(name)
	if err != nil {
   	 log.Fatal(err)
	}
	defer file.Close()
	/*scanner := bufio.NewScanner(file)
	
	for scanner.Scan() {
		str:=scanner.Text()
    	Dna=append(Dna,str)
	}
	
	if err := scanner.Err(); err != nil {
    	log.Fatal(err)
	}*/
	//r:= bufio.NewReaderSize(file, 4*1048576)
	// Some vectors are sequences of more than 30Kbp in one line...
	r:= bufio.NewReaderSize(file, 64536)
    
    line, isPrefix, err := r.ReadLine()
    for err == nil && !isPrefix {
        Dna=append(Dna,string(line))
        line, isPrefix, err = r.ReadLine()
    }
    if isPrefix {
        fmt.Println(errors.New("buffer size to small"))
        return Dna
    }
    if err != io.EOF {
        fmt.Println(err)
        return Dna
    }
	return Dna
}

func GetParametersFromline(line string) ([]int,error) {
	parm, err := ReadInts(strings.NewReader(line))
	return parm,err
}

func GetParametersFromline64(line string) ([]int64,error) {
	parm, err := ReadInts64(strings.NewReader(line))
	return parm,err
}