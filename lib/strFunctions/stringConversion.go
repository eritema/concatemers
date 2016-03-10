package strConversion

import (
   "strings"
   "strconv"
)


// Input: a string source that represents an array of strings combined with the separator string sep
// Output: an array of strings
func GetVectorOfStrings(source string, sep string) []string {
	return strings.Split(source,sep)
}

// Input: a string source that represents an array of int64 combined with the separator string sep
// Output: an array of int64
func GetVectorOfInt(source string, sep string) []int64 {
	
	vect:=strings.Split(source,sep)
	intVect:=make([]int64,len(vect))
	for i:=0;i<len(vect);i++ {
		intVect[i],_=strconv.ParseInt(vect[i],10,0)
	}
	return intVect
}	

// Input: a string s
// Output: the reverse of s
func ReverseString(s string) string {
    runes := []rune(s)
    for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
        runes[i], runes[j] = runes[j], runes[i]
    }
    return string(runes)
}