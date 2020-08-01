//
// Created by Craig Scott on 6/4/20.
//

#include <algorithm>
#include <vector>
#include <iostream>
#include <thread>

#include "vector.h"
#include "pcg32.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


constexpr float pi = 3.1415926535897932384626433832795f;
constexpr vec3f world_up = { 0, 1, 0 };
namespace intersect {
    constexpr float no_hit = -1;
}

struct ray
{
	vec3f o, d;
};


struct scene_object
{
	virtual vec3f get_color() const = 0;

	virtual ~scene_object() = default;

	virtual float intersect(const ray& r) const = 0;
	virtual vec3f get_normal(const vec3f& p_surface) const = 0;
};


struct sphere : public scene_object
{
    inline constexpr sphere() = default;
    inline constexpr sphere(float radius, vec3f origin, vec3f color) : m_radius(radius), m_origin(origin), m_color(color){}
    inline constexpr sphere(float radius) : m_radius(radius){}
    ~sphere() = default;

    inline const float get_radius() const {return m_radius;}
    inline const vec3f get_origin() const {return m_origin;}

    inline void set_radius(const float radius)  {m_radius = radius;}
    inline void set_origin(const float origin) {m_origin = origin;}

    vec3f get_color() const override{
        return m_color;
    }

    float intersect(const ray& r) const override{
        //Use geometric version of ray sphere intersection
        //This is the base of the right triangle with sides origin_ray-origin_sphere
        //ray_origin to sphere line bisect point
        //Define points Or (ray origin), Os(sphere origin).
        //Cast a ray into the sphere
        //Define points (A, B, C) for A first intersect, C for second intersect, B as the midpoint
        //Between A and C.
        //Define a line segment with starting point and end point in order AB, OrB, etc.
        //Iffy definitions but ok
        //

        //Find length from ray origin to sphere origin
        //Set ray origin is world origin

        const auto Or = r.o - m_origin;
        const auto Os = m_origin - m_origin;
        const auto OrOs = Os - Or;

        //Take dot product ray direction, to get component length, or length from Or to midpoint B
        const auto OrB_scalar = dot(OrOs,r.d);

        //Behind
        if(OrB_scalar < 0)
            return intersect::no_hit;

        const auto OsOr_scalar_squared = dot(OrOs,OrOs);
        const auto OrB_scalar_squared = OrB_scalar * OrB_scalar;
        const auto OsB_scalar_squared = OsOr_scalar_squared - OrB_scalar_squared;


        //Find length from bisect point to intersect point
        //Notice that OsA, OsB is just the radius of the sphere
        //We can calculate AB now or the point of intersection to midpoint
        const auto OsA_scalar = m_radius;
        const auto OsA_scalar_squared = OsA_scalar * OsA_scalar;
        const auto OsAOsB_squared_difference = OsA_scalar_squared - OsB_scalar_squared;

        //If short leg of right triangle is larger than radius, miss
        if(OsAOsB_squared_difference < 0)
            return intersect::no_hit;
        const auto AB_scalar = sqrt(OsAOsB_squared_difference);

        //return closest to ray origin
        return OrB_scalar - AB_scalar;


    }

    virtual vec3f get_normal(const vec3f& p_surface) const{
        vec3f surface_normal = p_surface - m_origin;
        surface_normal.normalize();
        return surface_normal;
    }
private:
    float m_radius{0};
    vec3f m_origin{0.0f,0.0f,0.0f};
    vec3f m_color{0.5f,0.5f,0.5f};
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
            if(t_int > 0 && t_int < t_min){
                t_min = t_int;
                obj_min = obj;
            }
        }
		return { obj_min,t_min };
	}
};


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


vec3f raytrace(const int pixel_x, const int pixel_y, const int sample_idx, int width, int height, const world const * world_ptr)
{
	const int pixel_idx = pixel_y * width + pixel_x;
	const vec2f pixel =
		vec2f((float)pixel_x,(float)pixel_y) +
		vec2f(randomized_halton(sample_idx, 0, pixel_idx),
		      randomized_halton(sample_idx, 1, pixel_idx));
	const float time = randomized_halton(sample_idx, 2, pixel_idx);

	const vec3f pinhole_pos{ 10.0f,10.0f,10.0f };
	const vec3f camera_dir = pinhole_pos * (-1.0f / pinhole_pos.magnitude());
	const vec3f camera_right = cross(camera_dir, world_up);
	const vec3f camera_up = cross(camera_right, camera_dir);
	const float fov = 80 * pi / 180;
	const float sensor_width = 2 * tanf(fov / 2);
	const float sensor_height = sensor_width * (float)height / width;
	const vec3f sensor_tlc = pinhole_pos + camera_dir - camera_right * (sensor_width / 2) + camera_up * (sensor_height / 2);
	const vec3f x_vec = camera_right * (sensor_width / width);
	const vec3f y_vec = camera_up * (-sensor_height / height);
	const vec3f sensor_pos = sensor_tlc + x_vec * pixel.x() + y_vec * pixel.y();
	const vec3f sensor_dir = sensor_pos - pinhole_pos;

	const ray r = { sensor_pos, sensor_dir * (1 / sensor_dir.magnitude()) };
    const auto scene_hit = world_ptr->nearest_intersection(r);
    const auto obj = scene_hit.first;
    const auto obj_t = scene_hit.second;
    if(obj) {
        const auto surface_pt = r.o + r.d * obj_t;
        const auto norm = obj->get_normal(surface_pt);
        const auto intensity = dot(r.d * -1, norm);
        return obj->get_color() * intensity;
    }
	return 0;
}


vec3f image_function(const int pixel_x, const int pixel_y, const int sample_idx, int width, int height, const world* const world_ptr)
{
#if 0
	const vec2f pixel = vec2f{ (float)pixel_x,(float)pixel_y } + vec2f{ halton(sample_idx, 0), halton(sample_idx,1) };
	const float time = halton(sample_idx, 2);
#else
	/*
	const int pixel_idx = pixel_y * width + pixel_x;
	const vec2f pixel =
		vec2f((float)pixel_x,(float)pixel_y) +
		vec2f(randomized_halton(sample_idx, 0 ,pixel_idx),
		      randomized_halton(sample_idx, 1, pixel_idx));
	//const float time = randomized_halton(sample_idx, 2, pixel_idx);
	 */
#endif
/*
	float circle = vec2f(pixel.x() - width / 2 + (time - 0.5f) * 240, pixel.y() - height / 2).magnitude();
	float r = circle > height / 4 ? 0.0f : 1.0f;
	float g = r;
	float b = r;
 */

	return raytrace(pixel_x, pixel_y, sample_idx, width, height, world_ptr);
}


int main(int argc, char* argv[])
{
	const int width = 512;
	const int height = 512;
	const int channels = 3;
	const int num_samples = 1 << 8;
    constexpr vec3f origin(0,0,0);
    constexpr vec3f default_color(0.5f,0.5f,0.5f);
	const std::unique_ptr<sphere> obj1 = std::make_unique<sphere>(3.0f, vec3f{1.0f,1.0f,1.0f}, default_color);

	//pcg32 rng;
    std::unique_ptr<world> world_container = std::make_unique<world>();
    world_container->objects.push_back(obj1.get());

	uint8_t* const pixels = new uint8_t[width * height * channels];
    std::atomic<int> curr_block(0);
	const int num_threads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads;
	threads.reserve(num_threads);
	constexpr int block_size = 16;
    const auto start_time = std::chrono::
	for(int thread_id = 0; thread_id < num_threads; ++thread_id) {
	    threads.push_back(std::thread([&]() {
	        //Initialize blocks
            while(true) {
                constexpr int x_blocks = (width + block_size - 1) / block_size; // Round up the number of blocks.
                constexpr int y_blocks = (height + block_size - 1) / block_size;

                const int block_idx = curr_block.fetch_add(1);

                if (block_idx >= x_blocks * y_blocks)
                    return;


                const int block_y  = block_idx / x_blocks;
                const int block_x  = block_idx - x_blocks * block_y; // Equiv to block_idx % x_blocks.
                const int block_x0 = block_x * block_size, block_x1 = std::min(block_x0 + block_size, width);
                const int block_y0 = block_y * block_size, block_y1 = std::min(block_y0 + block_size, height);

                for (int j = block_x0; j < block_x1; ++j) {
                    for (int i = block_y0; i < block_y1; ++i) {
                        vec3f sum = 0;
                        for (int s = 0; s < num_samples; ++s) {
                            const vec3f color_linear = raytrace(i, j, s, width, height, world_container.get());
                            sum += color_linear;
                        }
                        const vec3f color_linear = sum / num_samples;

                        const vec3f color_gamma(
                                powf(color_linear.x(), 1 / 2.2f),
                                powf(color_linear.y(), 1 / 2.2f),
                                powf(color_linear.z(), 1 / 2.2f));

                        const int ir = std::max(0, std::min(255, int(256 * color_gamma.x())));
                        const int ig = std::max(0, std::min(255, int(256 * color_gamma.y())));
                        const int ib = std::max(0, std::min(255, int(256 * color_gamma.z())));
                        const int index = width * j + i;
                        pixels[index * 3 + 0] = ir;
                        pixels[index * 3 + 1] = ig;
                        pixels[index * 3 + 2] = ib;
                    }
                }
            }
        }));
    }

	for(auto& t : threads){
	    t.join();
	}

	stbi_write_png("rt.png", width, height, channels, pixels, width * channels);

	delete [] pixels;

	return 0;
}
