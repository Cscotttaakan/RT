//
// Created by Craig Scott on 8/15/20.
//

#ifndef RT_UTILS_H
#define RT_UTILS_H


// Hash function by Thomas Wang: https://burtleburtle.net/bob/hash/integer.html
inline uint32_t hash(uint32_t x)
{
    x = (x ^ 12345391) * 2654435769;
    x ^= (x << 6) ^ (x >> 26);
    x *= 2654435769;
    x += (x << 5) ^ (x >> 12);

    return x;
}


// From PBRT
constexpr float FloatOneMinusEpsilon = 0.99999994f;
double RadicalInverse(int a, int base) noexcept
{
    const double invBase = 1.0 / base;

    int reversedDigits = 0;
    double invBaseN = 1;
    while (a)
    {
        const int next  = a / base;
        const int digit = a - base * next;
        reversedDigits  = reversedDigits * base + digit;
        invBaseN *= invBase;
        a = next;
    }
    return reversedDigits * invBaseN;
}


const static int max_primes = 6;
const static int primes[max_primes] = { 2, 3, 5, 7, 11, 13 };

inline float halton(int s, int d)
{
    return (float)std::min(RadicalInverse(s, primes[d % 6]), (double)FloatOneMinusEpsilon);
}

inline float randomized_halton(int s, int d, const int pixel_idx)
{
    const double u = RadicalInverse(s, primes[d % 6]);
    const double v = hash(pixel_idx * 17 + d) / 4294967296.0;
    const double w = (u + v < 1) ? u + v : u + v - 1;

    return (float)std::min(w, (double)FloatOneMinusEpsilon);
}
#endif //RT_UTILS_H
