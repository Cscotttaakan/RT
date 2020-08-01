//
// Created by Craig Scott on 6/6/20.
//

#pragma once

#include <array>
#include <cmath>



template<class TScalarType>
struct vector2
{
	TScalarType e[2];


	inline constexpr vector2() { }
	inline constexpr vector2(const vector2 & v) { for (int i = 0; i < 2; ++i) e[i] = v.e[i]; }
	inline constexpr vector2(const TScalarType & v) { for (int i = 0; i < 2; ++i) e[i] = v; }

	template<typename val, typename... vals, std::enable_if_t<(sizeof...(vals) > 0), int> = 0>
	inline constexpr vector2(const val v, const vals... vs) : e { (TScalarType)v, (TScalarType)vs... } { }

	inline constexpr TScalarType magnitude() const { return std::sqrt(e[0] * e[0] + e[1] * e[1]); }

	inline const TScalarType& x() const { return e[0]; }
	inline const TScalarType& y() const { return e[1]; }

	inline TScalarType& x() { return e[0]; }
	inline TScalarType& y() { return e[1]; }

	inline constexpr vector2 operator+(const vector2& rhs) const { return vector2(x() + rhs.x(), y() + rhs.y()); }
	inline constexpr vector2 operator-(const vector2& rhs) const { return vector2(x() - rhs.x(), y() - rhs.y()); }

	inline constexpr vector2& operator+=(const vector2& rhs) { e[0] += rhs.e[0]; e[1] += rhs.e[1]; return *this; }
	inline constexpr vector2& operator-=(const vector2& rhs) { e[0] -= rhs.e[0]; e[1] -= rhs.e[1]; return *this; }

	inline constexpr vector2 operator*(const TScalarType rhs) const { return vector3(x() * rhs, y() * rhs); }
	inline constexpr vector2 operator/(const TScalarType rhs) const { const TScalarType inv = 1 / rhs; return vector3(x() * inv, y() * inv); }
};


template <class TScalarType>
struct vector3
{
	TScalarType e[3];


	inline constexpr vector3() { }
	inline constexpr vector3(const vector3 & v) { for (int i = 0; i < 3; ++i) e[i] = v.e[i]; }
	inline constexpr vector3(const TScalarType & v) { for (int i = 0; i < 3; ++i) e[i] = v; }

	template<typename val, typename... vals, std::enable_if_t<(sizeof...(vals) > 0), int> = 0>
	inline constexpr vector3(const val v, const vals... vs) : e { (TScalarType)v, (TScalarType)vs... } { }

	inline constexpr TScalarType magnitude() const { return std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]); }
	inline constexpr void normalize() {const auto inverse_mag = 1/magnitude(); e[0]*=inverse_mag; e[1]*=inverse_mag; e[2]*=inverse_mag;}

	inline const TScalarType x() const { return e[0]; }
	inline const TScalarType y() const { return e[1]; }
	inline const TScalarType z() const { return e[2]; }

	inline TScalarType& x() { return e[0]; }
	inline TScalarType& y() { return e[1]; }
	inline TScalarType& z() { return e[2]; }

	inline constexpr vector3 operator-() const { return vector3(-x(), -y(), -z()); }

	inline constexpr vector3 operator+(const vector3& rhs) const { return vector3(x() + rhs.x(), y() + rhs.y(), z() + rhs.z()); }
	inline constexpr vector3 operator-(const vector3& rhs) const { return vector3(x() - rhs.x(), y() - rhs.y(), z() - rhs.z()); }

	inline constexpr vector3& operator+=(const vector3& rhs) { e[0] += rhs.e[0]; e[1] += rhs.e[1]; e[2] += rhs.e[2]; return *this; }
	inline constexpr vector3& operator-=(const vector3& rhs) { e[0] -= rhs.e[0]; e[1] -= rhs.e[1]; e[2] -= rhs.e[2]; return *this; }

	inline constexpr vector3 operator*(const TScalarType rhs) const { return vector3(x() * rhs, y() * rhs, z() * rhs); }
	inline constexpr vector3 operator/(const TScalarType rhs) const { const TScalarType inv = 1 / rhs; return vector3(x() * inv, y() * inv, z() * inv); }
};


template<class TScalarType>
TScalarType dot(const vector3<TScalarType>& a, const vector3<TScalarType>& b)
{
	return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}


template<class TScalarType>
vector3<TScalarType> cross(const vector3<TScalarType>& a, const vector3<TScalarType>& b)
{
	return vector3<TScalarType>(
		  a.y() * b.z() - b.y() * a.z(),
		-(a.x() * b.z() - b.x() * a.z()),
		  a.x() * b.y() - b.x() * a.y());
}


using vec2f = vector2<float>;
using vec3f = vector3<float>;

using vec2d = vector2<double>;
using vec3d = vector3<double>;
