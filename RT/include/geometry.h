//
// Created by Craig Scott on 8/15/20.
//

#ifndef RT_GEOMETRY_H
#define RT_GEOMETRY_H
#include "tiny_obj_loader.h"
#include "spectrum.h"
#include "vector.h"
struct ray
{
    vec3f o, d;
};


struct scene_object
{
    scene_object() = default;
    scene_object(const spectrum& albedo_, const spectrum& emission_) : albedo(albedo_), emission(emission_) { }
    virtual ~scene_object() = default;

    spectrum albedo = 0;
    spectrum emission = 0;
    float ior = 1.55f;
    bool mirror = false;

    virtual float intersect(const ray& r) const = 0;
    virtual vec3f get_normal(const vec3f& p_surface) const = 0;
};

using vert = vec3f;
using face = vec3i;
using normal = vec3f;

struct tri_mesh : public scene_object{
    inline tri_mesh() = default;
    virtual ~tri_mesh() = default;
    inline tri_mesh(const spectrum& albedo_, const spectrum& emission_, std::vector<face>&& faces):  scene_object(albedo_, emission_), m_faces(std::move(faces)){}
    float intersect(const ray& r){
        for(auto i = 0; i < m_faces.size(); ++i){
            //Backface culling
            if(dot(r.d, m_face_normals[i]) > 0.0f){

            }
        }
    }

    static bool triangle_ray_intersection(const ray& r, const vert& v1, const vert& v2, const vert& v3){
        //Treat v1 as origin for cross product
        //Vectors for all 3 sides
        const auto v1v2 = v2-v1;
        const auto v1v3 = v3-v1;
        const auto v3v2 = v2-v3;
        const float tri_area = (cross(v2-v1, v3-v1).magnitude())/2;
        const auto center = (v1 + v2 + v3)/3;


    }
private:
    std::vector<face> m_faces;
    std::vector<normal> m_face_normals;
    std::vector<vert> m_vertices;
};

struct sphere : public scene_object
{
    inline sphere() = default;
    inline sphere(const spectrum& albedo_, const spectrum& emission_, float radius, const vec3f& origin) : scene_object(albedo_, emission_), m_radius(radius), m_origin(origin) { }
    inline sphere(float radius) : m_radius(radius){}
    virtual ~sphere() = default;

    inline const float get_radius() const {return m_radius;}
    inline const vec3f get_origin() const {return m_origin;}

    inline void set_radius(const float radius)  {m_radius = radius;}
    inline void set_origin(const float origin) {m_origin = origin;}

    float intersect(const ray& r) const override
    {
        const vec3f s = r.o - m_origin;
        const float b = dot(s, r.d);
        const float c = dot(s, s) - m_radius * m_radius;

        const float discriminant = b * b - c;
        if (discriminant < 0)
            return -1;

        // Compute both roots and return the nearest one that's > 0
        const float t1 = -b - std::sqrt(discriminant);
        const float t2 = -b + std::sqrt(discriminant);
        return (t1 >= 3e-5f) ? t1 : t2;
    }

    virtual vec3f get_normal(const vec3f& p_surface) const override
    {
        return (p_surface - m_origin) / m_radius;
    }
private:
    float m_radius = 1;
    vec3f m_origin = 0;
};


struct world
{
    std::vector<scene_object*> objects;


    std::pair<scene_object*, float> nearest_intersection(const ray& r) const
    {
        float t_min = 10000000;
        scene_object* obj_min = nullptr;
        //Loop over all objects calling intersect, if distance > 0 && distance < t_min, update t_min and object min
        for(const auto& obj : objects){
            const auto t_int = obj->intersect(r);
            if(t_int > ray_eps && t_int < t_min){
                t_min = t_int;
                obj_min = obj;
            }
        }
        return { obj_min,t_min };
    }
};

#endif //RT_GEOMETRY_H
