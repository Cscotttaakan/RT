#### Notes

Loop over each pixel, check color. Can parallelize computations.

Looping in order + gamma correction.

#### Gamma

8 bits in RGB colors - Humans much more sensitive to darkness in image. Non-linearity in bits. Does not store
brightness linearly. 50% * 2 = 100%. Mathematically we work in linear space, when displaying we need logarithmic.

PNG loader -> bytes are going to be non-linear space (gamma space), we have to inverse gamma to get linear space.

Convert to linear space x^(2.2), convert from linear to sRGB x^(1/2.2)

Data is encoded in 8-bit numbers linearly within computation space, visually encoded logarithmically.

#### Monte-Carlo Integration

Randomly sample to find average height * base for integral area vs. classical integration methods
does not require as much sample points relative to dimensionality of integrand.


#### Quasi Monte-Carlo
Principle: For every sample s, there are d dimensions (x,y,t). For each we took a pure random function.
8 Binary Numbers / Reversed / Decimal Fraction of Second COL (1/2, 1/4, 1/8):
000 : 000 : 0.0
001 : 100 : 0.5
010 : 010 : 0.25
011 : 110 : 0.75
100 : 001 : 0.125
101 : 101 : 0.625
110 : 011 : 0.375
111 : 111 : 0.875

Can do same with base 3,5,... primes
Given index s, Halton Sequence defined as [RI_2(s), RI_3(s), RI_5(s),...,]

#### Ray-sphere intersection

Ray r(t) = o + d*t

Sphere at origin x^2 + y^2 + z^2 = r^2

x = o_x + d_x*t

x^2 = o_x^2 + 2o_x*d_x*t + (d_x*t)^2

o_x^2 + 2o_x*d_x*t + (d_x*t)^2 + o_y^2 + 2o_y*d_y*t + (d_y*t)^2 + o_z^2 + 2o_z*d_x*t + (d_z*t)^2 = r^2

(d_x^2 + d_y^2 + d_z^2)t^2 + 2(o_x*d_x + o_y*d_y + o_z*d_z)t + (o_x^2 + o_y^2 + o_z^2 - r^2) = 0

At^2 + Bt + C = 0
A = (d_x^2 + d_y^2 + d_z^2) = |d|^2 =  1 (Normalized vector)
B = 2 * dot(o,d)
C = dot(o,o) - r^2
B^2 - 4AC < 0 Ray misses sphere
(- B +- sqrt(B^2 - 4AC))/2A
Take negative sqrt discriminant to be closest


#### Questions - Readings

Montecarlo integration