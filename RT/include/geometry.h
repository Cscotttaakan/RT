//
// Created by Craig Scott on 8/15/20.
//

#ifndef RT_GEOMETRY_H
#define RT_GEOMETRY_H
#include <filesystem>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include "spectrum.h"
#include "vector.h"
constexpr float ray_eps = 5e-5f;
constexpr float T_MAX = 10000000;
constexpr int NUM_TRI_VERTS = 3;
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
using tex_coord = vec2f;
struct tri_mesh : public scene_object{
    inline tri_mesh() = default;
    virtual ~tri_mesh() = default;
    inline tri_mesh(const spectrum& albedo_, const spectrum& emission_,
            std::vector<face>&& faces,
            std::vector<vert>&& vertices,
            std::vector<normal>&& normals,
            std::vector<tex_coord>&& tex_coords)
            :  scene_object(albedo_, emission_),
            m_faces(std::move(faces)),
            m_vertices(std::move(vertices)),
            m_vert_normals(std::move(normals)),
            m_tex_coords(std::move(tex_coords)){}
    virtual float intersect(const ray& r) const override{
        float t_min = T_MAX;
        for(auto i = 0; i < m_faces.size(); ++i){
            auto t_val = T_MAX;
            auto intersect = triangle_ray_intersection(r,m_vertices[m_faces[i].x()], m_vertices[m_faces[i].y()], m_vertices[m_faces[i].z()], t_val);
            if(intersect && t_val - t_min < ray_eps)
                t_min = t_val;
        }
        return t_min;
    }

    virtual vec3f get_normal(const vec3f& p_surface) const override{
        return vec3f{};
    }


    static std::vector<tri_mesh> load_objs(const std::filesystem::path& obj_path){
        if(std::filesystem::exists(obj_path))
        {
            //Taken from example in Tinyobj repo
            tinyobj::attrib_t attrib;
            std::vector<tinyobj::shape_t> shapes;
            std::vector<tinyobj::material_t> materials;

            std::string warn;
            std::string err;

            //Default triangulate
            bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, obj_path.string().c_str());
            if (!warn.empty()) {
                std::cout << warn << std::endl;
            }

            if (!err.empty()) {
                std::cerr << err << std::endl;
            }
            std::vector<tri_mesh> meshes;
            meshes.reserve(shapes.size());
            // Loop over shapes
            for (size_t s = 0; s < shapes.size(); s++) {
                //Prepare data to be put into tri_mesh
                std::vector<face> faces;
                std::vector<vert> vertices;
                std::vector<normal> normals;
                std::vector<tex_coord> tex_coords;

                faces.reserve(shapes[s].mesh.num_face_vertices.size());
                vertices.reserve(NUM_TRI_VERTS * shapes[s].mesh.num_face_vertices.size());
                normals.reserve(NUM_TRI_VERTS * shapes[s].mesh.num_face_vertices.size());
                tex_coords.reserve(NUM_TRI_VERTS * shapes[s].mesh.num_face_vertices.size());
                size_t index_offset = 0;
                size_t vertex_index = 0;
                for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
                    int fv = shapes[s].mesh.num_face_vertices[f];
                    face face_value{};
                    // Loop over vertices in the face.
                    for (size_t v = 0; v < fv; v++) {
                        vert vertex_value{};
                        normal  normal_value{};
                        tex_coord tex_coord_value{};
                        //add face
                        face_value.e[v] = vertex_index;

                        // access to vertex
                        tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                        vertex_value.e[0] = attrib.vertices[3*idx.vertex_index+0];
                        vertex_value.e[1] = attrib.vertices[3*idx.vertex_index+1];
                        vertex_value.e[2] = attrib.vertices[3*idx.vertex_index+2];
                        normal_value.e[0] = attrib.normals[3*idx.normal_index+0];
                        normal_value.e[1] = attrib.normals[3*idx.normal_index+1];
                        normal_value.e[2] = attrib.normals[3*idx.normal_index+2];
                        tex_coord_value.e[0] = attrib.texcoords[2*idx.texcoord_index+0];
                        tex_coord_value.e[1] = attrib.texcoords[2*idx.texcoord_index+1];
                        // Optional: vertex colors
                        // tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
                        // tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
                        // tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];
                        vertices.push_back(std::move(vertex_value));
                        normals.push_back(std::move(normal_value));
                        tex_coords.push_back(std::move(tex_coord_value));
                        vertex_index++;
                    }
                    faces.push_back(std::move(face_value));
                    index_offset += fv;

                    // per-face material
                    shapes[s].mesh.material_ids[f];
                }
                tri_mesh mesh{1,
                              0,
                              std::move(faces),
                              std::move(vertices),
                              std::move(normals),
                              std::move(tex_coords)
                };
                meshes.push_back(std::move(mesh));
            }
            return meshes;
        }
        return std::vector<tri_mesh>{};
    }

    static bool triangle_ray_intersection(const ray& r, const vert& v1, const vert& v2, const vert& v3, float& tdist){
        // Center = u*v1 + v*v2 + w*
        const auto e1 = v2-v1;
        const auto e2 = v3-v1;
        const auto p = cross(r.d, e2);
        const float a = dot(e1, p);

        //parallel
        if(a < ray_eps && a > -ray_eps)
            return false;

        const float f = 1.0f/a;

        const vec3f s = r.o - v1;
        const float u = f * dot(s,p);
        if(u < 0.0f || u > 1.0f)
            return false;

        const vec3f q = dot(s,e1);
        const float v = f*dot(r.d,q);
        if(v < 0.0f || u + v > 1.0f)
            return false;
        auto tval = f*dot(e2,q);
        tdist = f*dot(e2,q);

        return tdist >= 0;
    }
private:
    std::vector<face> m_faces;
    std::vector<normal> m_vert_normals;
    std::vector<vert> m_vertices;
    std::vector<tex_coord> m_tex_coords;
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
        float t_min = T_MAX;
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
