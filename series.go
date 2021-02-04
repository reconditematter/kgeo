// Copyright (c) 2021 Leonid Kneller. All rights reserved.
// Licensed under the MIT license.
// See the LICENSE file for full license information.

package kgeo

// See https://geographiclib.sourceforge.io/html/geodseries30.html

func seriesA1(ε float64) float64 {
	const (
		c2 = 1.0 / 4.0
		c4 = 1.0 / 64.0
		c6 = 1.0 / 256.0
		c8 = 25.0 / 16384.0
	)
	ε2 := ε * ε
	return (1 + ε2*(c2+ε2*(c4+ε2*(c6+ε2*(c8))))) / (1 - ε)
}

func seriesA2(ε float64) float64 {
	const (
		c2 = 3.0 / 4.0
		c4 = 7.0 / 64.0
		c6 = 11.0 / 256.0
		c8 = 375.0 / 16384.0
	)
	ε2 := ε * ε
	return (1 - ε2*(c2+ε2*(c4+ε2*(c6+ε2*(c8))))) / (1 + ε)
}

func seriesA3(n, ε float64) float64 {
	c1 := 1.0/2.0 - 1.0/2.0*n
	c2 := 1.0/4.0 + n*(1.0/8.0+n*(-3.0/8.0))
	c3 := 1.0/16.0 + n*(3.0/16.0+n*(1.0/16.0+n*(-5.0/16.0)))
	c4 := 3.0/64.0 + n*(1.0/32.0+n*(5.0/32.0+n*(5.0/128.0+n*(-35.0/128.0))))
	c5 := 3.0/128.0 + n*(5.0/128.0+n*(5.0/256.0+n*(35.0/256.0+n*(7.0/256.0))))
	c6 := 5.0/256.0 + n*(15.0/1024.0+n*(35.0/1024.0+n*(7.0/512.0+n*(63.0/512.0))))
	c7 := 25.0/2048.0 + n*(35.0/2048.0+n*(21.0/2048.0+n*(63.0/2048.0+n*(21.0/2048.0))))
	c8 := 175.0/16384.0 + n*(35.0/4096.0+n*(63.0/4096.0+n*(63.0/8192.0+n*(231.0/8192.0))))
	return 1 - (ε * (c1 + ε*(c2+ε*(c3+ε*(c4+ε*(c5+ε*(c6+ε*(c7+ε*(c8)))))))))
}

func seriesC1(ε float64) (C1 [8]float64) {
	ε2 := ε * ε
	ε3 := ε * ε2
	ε4 := ε * ε3
	ε5 := ε * ε4
	ε6 := ε * ε5
	ε7 := ε * ε6
	ε8 := ε * ε7
	{
		const (
			c1 = -1.0 / 2.0
			c3 = +3.0 / 16.0
			c5 = -1.0 / 32.0
			c7 = +19.0 / 2048.0
		)
		C1[0] = ε * (c1 + ε2*(c3+ε2*(c5+ε2*(c7))))
	}
	{
		const (
			c2 = -1.0 / 16.0
			c4 = +1.0 / 32.0
			c6 = -9.0 / 2048.0
			c8 = +7.0 / 4096.0
		)
		C1[1] = ε2 * (c2 + ε2*(c4+ε2*(c6+ε2*(c8))))
	}
	{
		const (
			c3 = -1.0 / 48.0
			c5 = +3.0 / 256.0
			c7 = -3.0 / 2048.0
		)
		C1[2] = ε3 * (c3 + ε2*(c5+ε2*(c7)))
	}
	{
		const (
			c4 = -5.0 / 512.0
			c6 = +3.0 / 512.0
			c8 = -11.0 / 16384.0
		)
		C1[3] = ε4 * (c4 + ε2*(c6+ε2*(c8)))
	}
	{
		const (
			c5 = -7.0 / 1280.0
			c7 = +7.0 / 2048.0
		)
		C1[4] = ε5 * (c5 + ε2*(c7))
	}
	{
		const (
			c6 = -7.0 / 2048.0
			c8 = +9.0 / 4096.0
		)
		C1[5] = ε6 * (c6 + ε2*(c8))
	}
	{
		const (
			c7 = -33.0 / 14336.0
		)
		C1[6] = ε7 * (c7)
	}
	{
		const (
			c8 = -429.0 / 262144.0
		)
		C1[7] = ε8 * (c8)
	}
	return
}

func seriesC1p(ε float64) (C1p [8]float64) {
	ε2 := ε * ε
	ε3 := ε * ε2
	ε4 := ε * ε3
	ε5 := ε * ε4
	ε6 := ε * ε5
	ε7 := ε * ε6
	ε8 := ε * ε7
	{
		const (
			c1 = +1.0 / 2.0
			c3 = -9.0 / 32.0
			c5 = +205.0 / 1536.0
			c7 = -4879.0 / 73728.0
		)
		C1p[0] = ε * (c1 + ε2*(c3+ε2*(c5+ε2*(c7))))
	}
	{
		const (
			c2 = +5.0 / 16.0
			c4 = -37.0 / 96.0
			c6 = +1335.0 / 4096.0
			c8 = -86171.0 / 368640.0
		)
		C1p[1] = ε2 * (c2 + ε2*(c4+ε2*(c6+ε2*(c8))))
	}
	{
		const (
			c3 = +29.0 / 96.0
			c5 = -75.0 / 128.0
			c7 = +2901.0 / 4096.0
		)
		C1p[2] = ε3 * (c3 + ε2*(c5+ε2*(c7)))
	}
	{
		const (
			c4 = +539.0 / 1536.0
			c6 = -2391.0 / 2560.0
			c8 = +1082857.0 / 737280.0
		)
		C1p[3] = ε4 * (c4 + ε2*(c6+ε2*(c8)))
	}
	{
		const (
			c5 = +3467.0 / 7680.0
			c7 = -28223.0 / 18432.0
		)
		C1p[4] = ε5 * (c5 + ε2*(c7))
	}
	{
		const (
			c6 = +38081.0 / 61440.0
			c8 = -733437.0 / 286720.0
		)
		C1p[5] = ε6 * (c6 + ε2*(c8))
	}
	{
		const (
			c7 = +459485.0 / 516096.0
		)
		C1p[6] = ε7 * (c7)
	}
	{
		const (
			c8 = +109167851.0 / 82575360.0
		)
		C1p[7] = ε8 * (c8)
	}
	return
}

func seriesC2(ε float64) (C2 [8]float64) {
	ε2 := ε * ε
	ε3 := ε * ε2
	ε4 := ε * ε3
	ε5 := ε * ε4
	ε6 := ε * ε5
	ε7 := ε * ε6
	ε8 := ε * ε7
	{
		const (
			c1 = +1.0 / 2.0
			c3 = +1.0 / 16.0
			c5 = +1.0 / 32.0
			c7 = +41.0 / 2048.0
		)
		C2[0] = ε * (c1 + ε2*(c3+ε2*(c5+ε2*(c7))))
	}
	{
		const (
			c2 = +3.0 / 16.0
			c4 = +1.0 / 32.0
			c6 = +35.0 / 2048.0
			c8 = +47.0 / 4096.0
		)
		C2[1] = ε2 * (c2 + ε2*(c4+ε2*(c6+ε2*(c8))))
	}
	{
		const (
			c3 = +5.0 / 48.0
			c5 = +5.0 / 256.0
			c7 = +23.0 / 2048.0
		)
		C2[2] = ε3 * (c3 + ε2*(c5+ε2*(c7)))
	}
	{
		const (
			c4 = +35.0 / 512.0
			c6 = +7.0 / 512.0
			c8 = +133.0 / 16384.0
		)
		C2[3] = ε4 * (c4 + ε2*(c6+ε2*(c8)))
	}
	{
		const (
			c5 = +63.0 / 1280.0
			c7 = +21.0 / 2048.0
		)
		C2[4] = ε5 * (c5 + ε2*(c7))
	}
	{
		const (
			c6 = +77.0 / 2048.0
			c8 = +33.0 / 4096.0
		)
		C2[5] = ε6 * (c6 + ε2*(c8))
	}
	{
		const (
			c7 = +429.0 / 14336.0
		)
		C2[6] = ε7 * (c7)
	}
	{
		const (
			c8 = +6435.0 / 262144.0
		)
		C2[7] = ε8 * (c8)
	}
	return
}
