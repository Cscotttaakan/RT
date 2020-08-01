#pragma once


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
};
