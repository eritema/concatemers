package lib

import (
    "bufio"
    "io"
    "strconv"
)

// ReadInts reads whitespace-separated ints from r. If there's an error, it
// returns the ints successfully read so far as well as the error value.
// Example
// ReadInts(strings.NewReader(string))
func ReadInts(r io.Reader) ([]int, error) {
    scanner := bufio.NewScanner(r)
    scanner.Split(bufio.ScanWords)
    var result []int
    var x int
    var err error
    for scanner.Scan() {
        x, err = strconv.Atoi(scanner.Text())
        if err != nil {
            return result, err
        }
        result = append(result, x)
    }
    return result, scanner.Err()
}

func ReadInts64(r io.Reader) ([]int64, error) {
    scanner := bufio.NewScanner(r)
    scanner.Split(bufio.ScanWords)
    var result []int64
    var x int64
    var err error
    for scanner.Scan() {
        x, err = strconv.ParseInt(scanner.Text(),10,0)
        if err != nil {
            return result, err
        }
        result = append(result, x)
    }
    return result, scanner.Err()
}



func Readfloats64(r io.Reader) ([]float64, error) {
    scanner := bufio.NewScanner(r)
    scanner.Split(bufio.ScanWords)
    var result []float64
    for scanner.Scan() {
        x, err := strconv.ParseFloat(scanner.Text(),64)
        if err != nil {
            return result, err
        }
        result = append(result, x)
    }
    return result, scanner.Err()
}