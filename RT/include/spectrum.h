#pragma once
#include "vector.h"

constexpr double wavelength_min = 380.0f;
constexpr double wavelength_max = 740.0f;

struct spectrum
{
	float e[3];

	inline spectrum() { }
	inline spectrum(const spectrum & v) { for (int i = 0; i < 3; ++i) e[i] = v.e[i]; }
	inline spectrum(const float & v) { for (int i = 0; i < 3; ++i) e[i] = v; }

	inline spectrum(float x, float y, float z) { e[0] = x; e[1] = y; e[2] = z; }

	//template<typename val, typename... vals, std::enable_if_t<(sizeof...(vals) > 0), int> = 0>
	//inline constexpr spectrum(const val v, const vals... vs) : e { (float)v, (float)vs... } { }

	inline const float x() const { return e[0]; }
	inline const float y() const { return e[1]; }
	inline const float z() const { return e[2]; }

	inline float& x() { return e[0]; }
	inline float& y() { return e[1]; }
	inline float& z() { return e[2]; }

	inline spectrum operator+(const spectrum& rhs) const { return spectrum(x() + rhs.x(), y() + rhs.y(), z() + rhs.z()); }
	inline spectrum operator-(const spectrum& rhs) const { return spectrum(x() - rhs.x(), y() - rhs.y(), z() - rhs.z()); }

	inline spectrum& operator+=(const spectrum& rhs) { e[0] += rhs.e[0]; e[1] += rhs.e[1]; e[2] += rhs.e[2]; return *this; }
	inline spectrum& operator-=(const spectrum& rhs) { e[0] -= rhs.e[0]; e[1] -= rhs.e[1]; e[2] -= rhs.e[2]; return *this; }

	inline spectrum operator*(const float rhs) const { return spectrum(x() * rhs, y() * rhs, z() * rhs); }
	inline spectrum operator/(const float rhs) const { const float inv = 1 / rhs; return spectrum(x() * inv, y() * inv, z() * inv); }

	inline spectrum operator *(const spectrum & rhs) const { spectrum s; for (int i = 0; i < 3; ++i) s.e[i] = e[i] * rhs.e[i]; return s; }

	inline spectrum& operator*=(const spectrum & rhs) { for (int i = 0; i < 3; ++i) e[i] *= rhs.e[i]; return *this; }

	static vec3f wavelength2xyz(float wave){
	    vec3f spec = 0;
	    spec.x() = xFit_1931(wave);
	    spec.y() = yFit_1931(wave);
	    spec.z() = zFit_1931(wave);
	    return spec;
	}

	static vec3f xyz2rgb(const vec3f& xyz){
	    vec3f color = 0;
        color.x() =  3.2404542*xyz.x() - 1.5371385*xyz.y() - 0.4985314*xyz.z();
        color.y() = -0.9692660*xyz.x() + 1.8760108*xyz.y() + 0.0415560*xyz.z();
        color.z() =  0.0556434*xyz.x() - 0.2040259*xyz.y() + 1.0572252*xyz.z();
        return std::move(color);
	}

	static vec3f rgb2wavelength(const vec3f& rgb){

	}

    static float xFit_1931( float wave )
    {
        float t1 = (wave-442.0f)*((wave<442.0f)?0.0624f:0.0374f);
        float t2 = (wave-599.8f)*((wave<599.8f)?0.0264f:0.0323f);
        float t3 = (wave-501.1f)*((wave<501.1f)?0.0490f:0.0382f);
        return 0.362f*expf(-0.5f*t1*t1) + 1.056f*expf(-0.5f*t2*t2)
               - 0.065f*expf(-0.5f*t3*t3);
    }
    static float yFit_1931( float wave )
    {
        float t1 = (wave-568.8f)*((wave<568.8f)?0.0213f:0.0247f);
        float t2 = (wave-530.9f)*((wave<530.9f)?0.0613f:0.0322f);
        return 0.821f*exp(-0.5f*t1*t1) + 0.286f*expf(-0.5f*t2*t2);
    }
    static float zFit_1931( float wave )
    {
        float t1 = (wave-437.0f)*((wave<437.0f)?0.0845f:0.0278f);
        float t2 = (wave-459.0f)*((wave<459.0f)?0.0385f:0.0725f);
        return 1.217f*exp(-0.5f*t1*t1) + 0.681f*expf(-0.5f*t2*t2);
    }

};
