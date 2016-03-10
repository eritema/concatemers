package lib

import (	
	//"fmt"
)

// Cartesian Point

type Point struct{
	X,Y int
}

// vector of 2 components
type Vector struct {
	X,Y,length int
}


//Rectangle defined with 3 points A,B,C
// AB and BC are the orthogonal edges
type Rect struct {
	A,B,C Point
}


//algebric dot product
func Dot(v1, v2 Vector) int {
	return v1.X*v2.X+v1.Y*v2.Y
}


// Given 2 points A and B return the vector AB
func Vectorize (p1,p2 Point) Vector {
	v:=Vector{X:p2.X-p1.X, Y:p2.Y-p1.Y}
	v.length=Dot(v,v)
	return v
}

type RectSet struct {
	PointLists []Point
}


// Given a point p1, built two oblique rects (+-45 deg) around it and check if p2 is contained
// return an int flag: 0 -> no; 1 -> yes (+45); 2 -> yes (-45); 3 -> yes (+-45)


func InRectPlus(p1,p2 Point, w,d int) int {
	var R Rect
	var inFlag int
	inFlag=0
	R.A=Point{X:p1.X+d, Y:p1.Y-d}
	R.B=Point{X:p1.X-d, Y:p1.Y+d}
	R.C=Point{X:p1.X+w-d, Y:p1.Y+w+d}
	AB:=Vectorize(R.A,R.B)
	BC:=Vectorize(R.B,R.C)
	dotAB_AP:=Dot(AB,Vectorize(R.A,p2))
	dotBC_BP:=Dot(BC,Vectorize(R.B,p2))
	if AB.length>=dotAB_AP && dotAB_AP>=0 && BC.length>=dotBC_BP && dotBC_BP >= 0 {
		inFlag=1
	}
	return inFlag
}

func IsPointIn(point Point, pointSet []RectSet, strand string, w,d int) []int {	
	var fragments []int
	if strand!="minus" {
		for i:=0;i<len(pointSet);i++ {
			p:=pointSet[i]
			flag:=InRectPlus(p.PointLists[len(p.PointLists)-1],point, w,d)
			if flag==1 {
				fragments=append(fragments,i)
			} 	
		}
	} else {	
		for i:=0;i<len(pointSet);i++ {
			p:=pointSet[i]
			flag:=InRectMinus(p.PointLists[len(p.PointLists)-1],point, w,d)
			if flag==1 {
				fragments=append(fragments,i)
			} 
		}
	}
	return fragments
}



func InRectMinus(p1,p2 Point, w,d int) int {
	var Rc Rect
	var inFlag int
	inFlag=0
	Rc.A=Point{X:p1.X-d, Y:p1.Y-d}
	Rc.B=Point{X:p1.X+d, Y:p1.Y+d}
	Rc.C=Point{X:p1.X+w+d, Y:p1.Y-w+d}
	ABc:=Vectorize(Rc.A,Rc.B)
	BCc:=Vectorize(Rc.B,Rc.C)
	dotAB_APc:=Dot(ABc,Vectorize(Rc.A,p2))
	dotBC_BPc:=Dot(BCc,Vectorize(Rc.B,p2))
	if ABc.length>=dotAB_APc && dotAB_APc >=0 && BCc.length>=dotBC_BPc && dotBC_BPc >=0 {
		inFlag=1
	}
	return inFlag
}


