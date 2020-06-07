//
// Created by Craig Scott on 6/6/20.
//

#ifndef RT_VECTOR_H
#define RT_VECTOR_H
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

    virtual ~vector2() = default;

    TScalarType magnitude(){return sqrt(pow(m_data[0], 2) + pow(m_data[1], 2));}

    TScalarType x() { return m_data[0]; }

    TScalarType y() { return m_data[1]; }

    vector2<TScalarType> operator*(const vector2<TScalarType> &a) {
        return vector2(m_data[0] + a.x(), m_data[1] + a.y());
    }

    TScalarType* data() { return m_data.data(); };
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

    virtual ~vector3() = default;

    TScalarType magnitude(){return sqrt(pow(m_data[0], 2) + pow(m_data[1], 2) + pow(m_data[2], 2));}

    TScalarType x() { return m_data[0]; }

    TScalarType y() { return m_data[1]; }

    TScalarType z() { return m_data[2]; }

    vector2<TScalarType> operator*(const vector2<TScalarType> &a) {
        return vector2(m_data[0] + a.x(), m_data[1] + a.y(), m_data[2] + a.z());
    }

    vector2<TScalarType> cross(const vector2<TScalarType> &a) {
        return vector2(y()*a.z() - a.y()*z(), -(x() * a.z() - a.x()*z()), x() * a.y() - a.x()*y());
    }

    TScalarType* data() { return m_data.data(); };
private:
    std::array<TScalarType, 3> m_data;
};

using vec2f = vector2<float>;
using vec2d = vector2<double>;
using vec3f = vector2<float>;
using vec3d = vector2<double>;
#endif //RT_VECTOR_H
