package reconstruction


func Overlap(km1,km2 string) bool {
	if len(km1)!=len(km2) {
		return false
	}
	return km1[1:len(km1)]==km2[0:len(km1)-1]
}
