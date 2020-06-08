//
// Created by Craig Scott on 6/4/20.
//

#include <algorithm>
#include <vector>

#include "vector.h"
#include "stb_write_impl.h"
#include "pcg32.h"

using byte = unsigned char;
// From PBRT
constexpr double DoubleOneMinusEpsilon = 0.99999999999999989;
constexpr float  FloatOneMinusEpsilon  = 0.99999994f;
const vec3f world_up{0.0f, 1.0f,0.0};

struct ray{
    vec3f o,d;
};

struct scene_object{
    virtual ~scene_object() = default;

    virtual float intersect(const ray& r)const = 0;

    virtual vec3f get_normal(const vec3f& p_surface)const = 0;

    //virtual vec3f color
};

struct sphere : public scene_object{

};

struct world{
    std::vector<scene_object*> objects;

    std::pair<scene_object*,float> nearest_intersection(const ray& r) const{
        float t_min = 10000000;
        scene_object* obj_min = nullptr;
        //Loop over all objects calling intersect, if distance > 0 && distance < t_min, update t_min and object min

        return {obj_min,t_min};
    }

};

// Hash function by Thomas Wang: https://burtleburtle.net/bob/hash/integer.html
inline uint32_t hash(uint32_t x)
{
    x  = (x ^ 12345391) * 2654435769;
    x ^= (x << 6) ^ (x >> 26);
    x *= 2654435769;
    x += (x << 5) ^ (x >> 12);

    return x;
}

float RadicalInverse(int a, int base) noexcept
{
    const double invBase = 1.0 / base;

    int reversedDigits = 0;
    double invBaseN = 1;
    while (a)
    {
        const int next  = a / base;
        const int digit = a - base * next;
        reversedDigits = reversedDigits * base + digit;
        invBaseN *= invBase;
        a = next;
    }

    return (float)std::min(reversedDigits * invBaseN, (double)FloatOneMinusEpsilon);
}

float halton(int s, int d){
    const static int max_primes = 6;
    const static int primes[max_primes] = {2,3,5,7,11,13};
    return RadicalInverse(s,primes[d%6]);
}

float randomized_halton(int s, int d, const int pixel_idx){
    const static int max_primes = 6;
    const static int primes[max_primes] = {2,3,5,7,11,13};
    const float u = RadicalInverse(s,primes[d%6]);
    const float v = hash(pixel_idx*17 + d)/4294967296.0f;
    const float w = (u+v < 1) ? u+v : u+v - 1;
    return w;
}

vec3f raytrace(const int pixel_x, const int pixel_y, const int sample_idx, int width, int height){
    const int pixel_idx = pixel_y * width + pixel_x;
    const vec2f pixel = vec2f{(float)pixel_x,(float)pixel_y} +
                        vec2f{randomized_halton(sample_idx, 0,pixel_idx),
                              randomized_halton(sample_idx,1,pixel_idx)};
    const float time = randomized_halton(sample_idx, 2, pixel_idx);

    const vec3f pinhole_pos{2.0f,4.0f,1.0f};
    const vec3f camera_dir = pinhole_pos * (-1.0f/pinhole_pos.magnitude()) ;
    const vec3f camera_right = cross(camera_dir, world_up);
    const vec3f camera_up = cross(camera_right, camera_dir);
    const float fov = 80*M_PI/180;
    const float sensor_width = 2 * tanf(fov/2);
    const float sensor_height = sensor_width * (float)height/width;
    const vec3f sensor_tlc = pinhole_pos + camera_dir - camera_right * (sensor_width/2) + camera_up *(sensor_height/2);
    const vec3f x_vec = camera_right * (sensor_width/width);
    const vec3f y_vec = camera_up * (-sensor_height/height);
    const vec3f sensor_pos = sensor_tlc + x_vec*pixel.x() + y_vec * pixel.y();
    const vec3f sensor_dir = pinhole_pos - sensor_pos;

    const ray r = {sensor_pos, sensor_dir * (1/sensor_dir.magnitude())};
}

vec3f image_function(const int pixel_x, const int pixel_y, const int sample_idx, int width, int height){
#if 0
    const vec2f pixel = vec2f{(float)pixel_x,(float)pixel_y} + vec2f{halton(sample_idx, 0), halton(sample_idx,1)};
    const float time = halton(sample_idx, 2);
#else
    const int pixel_idx = pixel_y * width + pixel_x;
    const vec2f pixel = vec2f{(float)pixel_x,(float)pixel_y} +
            vec2f{randomized_halton(sample_idx, 0,pixel_idx),
            randomized_halton(sample_idx,1,pixel_idx)};
    const float time = randomized_halton(sample_idx, 2, pixel_idx);
#endif
    float circle = vec2f(pixel.x() - width/2 + (time - 0.5f)*240, pixel.y() - height/2).magnitude();
    float r = circle > height/4 ? 0.0f : 1.0f ;
    float g = r;
    float b = r;
    return vec3f(r,g,b);
}

int main(int argc, char *argv[])
{
    const int width = 512;
    const int height = 512;
    const int channels = 3;
    const int num_samples = 1 << 10;
    pcg32 rng;
    uint8_t* pixels = new uint8_t[width*height*channels];

#pragma omp parallel for
    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            vec3f sum{0,0,0};
            for(int s = 0; s < num_samples; ++s){
                const vec3f color_linear = image_function(i,j,s, width, height);
                sum += color_linear;
            }
            const vec3f color_linear = sum*(1.0f/num_samples);

            const vec3f color_gamma{
                powf(color_linear.x(),1/2.2f),
                powf(color_linear.y(),1/2.2f),
                powf(color_linear.z(),1/2.2f)};

            const int ir = std::max(0,std::min(255,int(255.99 * color_gamma.x())));
            const int ig = std::max(0,std::min(255,int(255.99 * color_gamma.y())));
            const int ib = std::max(0,std::min(255,int(255.99 * color_gamma.z())));
            const int index = width*j + i;
            pixels[index*3 + 0] = ir;
            pixels[index*3 + 1] = ig;
            pixels[index*3 + 2] = ib;
        }
    }

    stbi_write_png("64_sample_rqmc_motion_blur.png", width, height, channels, pixels, width * channels);
    delete [] pixels;

    return 0;
}
