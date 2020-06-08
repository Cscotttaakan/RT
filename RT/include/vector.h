//
// Created by Craig Scott on 6/6/20.
//

#pragma once

#include <array>
#include <cmath>

template<class TScalarType>
class vector2 {
public:
    vector2() = default;
    vector2(TScalarType x, TScalarType y) {
        m_data[0] = x;
        m_data[1] = y;
    }

    ~vector2() = default;

    const TScalarType magnitude() const {return sqrt(m_data[0]*m_data[0] + m_data[1]*m_data[1]);}

    const TScalarType x() const{ return m_data[0]; }

    const TScalarType y() const{ return m_data[1]; }

    TScalarType& x() { return m_data[0]; }

    TScalarType& y() { return m_data[1]; }
    TScalarType* data() { return m_data.data(); };

    vector2 operator +(const vector2& rhs) const{
        return vector2(x() + rhs.x(), y() + rhs.y());
    }
private:
    std::array<TScalarType, 2> m_data;
};

template <class TScalarType>
class vector3 {
public:
    vector3() = default;
    vector3(TScalarType x, TScalarType y, TScalarType z) {
        m_data[0] = x;
        m_data[1] = y;
        m_data[2] = z;
    }

    ~vector3() = default;

    const TScalarType magnitude() const{return std::sqrt(m_data[0]*m_data[0] + m_data[1]*m_data[1] + m_data[2]*m_data[2]);}

    const TScalarType x() const { return m_data[0]; }

    const TScalarType y() const { return m_data[1]; }

    const TScalarType z() const { return m_data[2]; }

    TScalarType& x() { return m_data[0]; }

    TScalarType& y() { return m_data[1]; }

    TScalarType& z() { return m_data[2]; }

    TScalarType* data() { return m_data.data(); };

    const vector3& operator+=(const vector3& rhs){
        m_data[0] += rhs.x();
        m_data[1] += rhs.y();
        m_data[2] += rhs.z();
        return *this;
    }

    const vector3 operator *(const TScalarType rhs) const{
        return vector3(rhs*x(), rhs*y(), rhs*z());
    }

    const vector3 operator +(const vector3& rhs) const{
        return vector3(rhs.x() + x(), rhs.y() + y(), rhs.z() + z());
    }

    const vector3 operator -(const vector3& rhs) const{
        return vector3(x() - rhs.x(), y() - rhs.y(), z() - rhs.z());
    }

private:
    std::array<TScalarType, 3> m_data;
};

template<class TScalarType>
TScalarType dot(const vector3<TScalarType> &a, const vector3<TScalarType> &b) {
    return b.x()*a.x() + b.y()*a.y() + b.z()*a.z();
}

template<class TScalarType>
vector3<TScalarType> cross(const vector3<TScalarType> &a, const vector3<TScalarType> &b) {
    return vector3<TScalarType>(b.y()*a.z() - a.y()*b.z(), -(b.x() * a.z() - a.x()*b.z()), b.x() * a.y() - a.x()*b.y());
}

/*
vector2<TScalarType> dot(const vector2<TScalarType> &a,const vector2<TScalarType> &b) {
    return vector2(b.x() + a.x(), b.x() + a.y());
}
*/
using vec2f = vector2<float>;
using vec2d = vector2<double>;
using vec3f = vector3<float>;
using vec3d = vector3<double>;
