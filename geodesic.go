// Copyright (c) 2021 Leonid Kneller. All rights reserved.
// Licensed under the MIT license.
// See the LICENSE file for full license information.

package kgeo

import (
	"math"
)

// Geodesic -- geodesic solver for a spheroidal model of the Earth.
//
// Reference: Karney, C.F.F. Algorithms for geodesics. J Geod 87, 43–55 (2013).
//
// DOI: https://doi.org/10.1007/s00190-012-0578-z
type Geodesic struct {
	a, b    float64
	f, n    float64
	e2, ep2 float64
	cA3     *coeffA3
	cC3     *coeffC3
}

// Solution -- represents a solution of direct/inverse geodesic problems.
type Solution struct {
	Lat1, Lon1, Azi1 float64
	Lat2, Lon2, Azi2 float64
	S12              float64
}

// NewGeodesic -- returns a geodesic solver for the spheroid
// defined by `a` (equatorial axis) and `f` (flattening).
func NewGeodesic(a, f float64) Geodesic {
	if !(1 <= a && a <= 1e10) {
		panic("kgeo.NewGeodesic: invalid argument `a`")
	}
	if !(0 <= f && f <= 1.0/150) {
		panic("kgeo.NewGeodesic: invalid argument `f`")
	}
	if f <= 1.0/(1<<26) {
		f = 0
	}
	g := Geodesic{a: a, b: a * (1 - f), f: f, n: f / (2 - f), e2: f * (2 - f), ep2: f * (2 - f) / ((1 - f) * (1 - f))}
	g.cA3 = newCoeffA3(g.n)
	g.cC3 = newCoeffC3(g.n)
	return g
}

// A -- returns the equatorial (major) axis (a) of `g`.
func (g Geodesic) A() float64 {
	return g.a
}

// F -- returns the (first) flattening (f) of `g`.
func (g Geodesic) F() float64 {
	return g.f
}

// Direct -- solves the direct problem: given the source point defined by `lat1` and `lon1`,
// the azimuth `azi1`, and the geodesic length `s12`, find the target point and the azimuth
// at that point.
func (g Geodesic) Direct(lat1, lon1, azi1, s12 float64) Solution {
	if !(-90 <= lat1 && lat1 <= +90) {
		panic("kgeo.Geodesic.Direct: invalid argument `lat1`")
	}
	if !(-180 <= lon1 && lon1 <= +180) {
		panic("kgeo.Geodesic.Direct: invalid argument `lon1`")
	}
	if !(-180 <= azi1 && azi1 <= +180) {
		panic("kgeo.Geodesic.Direct: invalid argument `azi1`")
	}
	if !(0 <= s12 && s12 <= 1e11) {
		panic("kgeo.Geodesic.Direct: invalid argument `s12`")
	}
	//
	{
		// adjust near-polar latitudes
		const ε = 1.0 / (1 << 38)
		if math.Abs(lat1) > 90*(1-ε) {
			lat1 = math.Copysign(90*(1-ε), lat1)
		}
	}
	//
	b := g.b
	f := g.f
	ep2 := g.ep2
	//
	φ1 := lat1 * (math.Pi / 180)
	sinφ1, cosφ1 := math.Sincos(φ1)
	α1 := azi1 * (math.Pi / 180)
	sinα1, cosα1 := math.Sincos(α1)
	// solve triangle NEA
	β1 := math.Atan2((1-f)*sinφ1, cosφ1)
	sinβ1, cosβ1 := math.Sincos(β1)
	α0 := math.Atan2(sinα1*cosβ1, math.Hypot(cosα1, sinα1*sinβ1))
	sinα0, cosα0 := math.Sincos(α0)
	σ1 := math.Atan2(sinβ1, cosα1*cosβ1)
	sinσ1, cosσ1 := math.Sincos(σ1)
	ω1 := math.Atan2(sinα0*sinσ1, cosσ1)
	// determine σ2
	var tt float64
	k2 := ep2 * cosα0 * cosα0
	tt = math.Sqrt(1 + k2)
	ε := (tt - 1) / (tt + 1)
	A1 := seriesA1(ε)
	C1 := seriesC1(ε)
	I1σ1 := A1 * (σ1 + sumSin(σ1, C1))
	s1 := b * I1σ1
	s2 := s1 + s12
	τ2 := s2 / (b * A1)
	C1p := seriesC1p(ε)
	σ2 := τ2 + sumSin(τ2, C1p)
	sinσ2, cosσ2 := math.Sincos(σ2)
	// solve triangle NEB
	α2 := math.Atan2(sinα0, cosα0*cosσ2)
	β2 := math.Atan2(cosα0*sinσ2, math.Hypot(cosα0*cosσ2, sinα0))
	sinβ2, cosβ2 := math.Sincos(β2)
	ω2 := math.Atan2(sinα0*sinσ2, cosσ2)
	φ2 := math.Atan2(sinβ2, (1-f)*cosβ2)
	// determine λ12
	A3 := g.cA3.seriesA3(ε)
	C3 := g.cC3.seriesC3(ε)
	I3σ1 := A3 * (σ1 + sumSin(σ1, C3))
	I3σ2 := A3 * (σ2 + sumSin(σ2, C3))
	λ1 := ω1 - f*sinα0*I3σ1
	λ2 := ω2 - f*sinα0*I3σ2
	λ12 := λ2 - λ1
	//
	lat2 := φ2 * (180 / math.Pi)
	lon2 := lon1 + λ12*(180/math.Pi)
	if lon2 < -180 {
		lon2 += 360
	}
	if lon2 > +180 {
		lon2 -= 360
	}
	azi2 := α2 * (180 / math.Pi)
	//
	return Solution{Lat1: nnz(lat1), Lon1: nnz(lon1), Azi1: nnz(azi1), Lat2: nnz(lat2), Lon2: nnz(lon2), Azi2: nnz(azi2), S12: nnz(s12)}
}

func (g Geodesic) hybrid(sinβ1, cosβ1, sinβ2, cosβ2, sinα1, cosα1 float64) (float64, float64) {
	b := g.b
	ep2 := g.ep2
	// solve triangle NEA
	α0 := math.Atan2(sinα1*cosβ1, math.Hypot(cosα1, sinα1*sinβ1))
	sinα0, cosα0 := math.Sincos(α0)
	σ1 := math.Atan2(sinβ1, cosα1*cosβ1)
	// solve triangle NEB
	α2 := math.Atan2(sinα0, math.Sqrt(sq(cosα1*cosβ1)+(cosβ2-cosβ1)*(cosβ2+cosβ1)))
	cosα2 := math.Cos(α2)
	σ2 := math.Atan2(sinβ2, cosα2*cosβ2)
	// determine s12 and λ12
	k2 := ep2 * cosα0 * cosα0
	tt := math.Sqrt(1 + k2)
	ε := (tt - 1) / (tt + 1)
	A1 := seriesA1(ε)
	C1 := seriesC1(ε)
	I1σ1 := A1 * (σ1 + sumSin(σ1, C1))
	s1 := b * I1σ1
	I1σ2 := A1 * (σ2 + sumSin(σ2, C1))
	s2 := b * I1σ2
	s12 := s2 - s1
	//
	return α2, s12
}

// sq -- square
func sq(x float64) float64 {
	return x * x
}

// nnz -- no negative zero
func nnz(x float64) float64 {
	if x == 0 {
		return 0
	}
	return x
}
