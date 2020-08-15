//
// Created by Craig Scott on 6/4/20.
//

#include <algorithm>
#include <atomic>
#include <vector>
#include <iostream>
#include <thread>

#include "vector.h"
#include "spectrum.h"

#include "pcg32.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


constexpr float float_pi = 3.1415926535897932384626433832795f;
constexpr float float_twopi = 6.283185307179586476925286766559f;
constexpr float ray_eps = 5e-5f;
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


vec3f raytrace(const int pixel_x, const int pixel_y, const int sample_idx, int width, int height, const world * const world_ptr)
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
	const float fov = 80 * float_pi / 180;
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
        //return obj->albedo * intensity;
        return intensity;
    }
	return 0;
}

double chooseRandWavelength(){
    double f = (double)rand() / RAND_MAX;
    return wavelength_min + f * (wavelength_max - wavelength_min);
}


vec3f pathtrace(const int pixel_x, const int pixel_y, const int sample_idx, int width, int height, const world * const world_ptr, pcg32 * rng)
{
    constexpr int max_bounces = 32;

    const int pixel_idx = pixel_y * width + pixel_x;
    const vec2f pixel =
        vec2f((float)pixel_x,(float)pixel_y) +
        vec2f(
            randomized_halton(sample_idx, 0, pixel_idx),
            randomized_halton(sample_idx, 1, pixel_idx));
    const vec3f pinhole_pos{ 10.0f,8.0f,10.0f };
    vec3f camera_dir = -pinhole_pos + vec3f{ 0, 3, 0};
    camera_dir.normalize();
    const vec3f camera_right = cross(camera_dir, world_up);
    const vec3f camera_up = cross(camera_right, camera_dir);
    const float fov = 30 * float_pi / 180;
    const float sensor_width = 2 * tanf(fov / 2);
    const float sensor_height = sensor_width * (float)height / width;
    const vec3f sensor_tlc = pinhole_pos + camera_dir - camera_right * (sensor_width / 2) + camera_up * (sensor_height / 2);
    const vec3f x_vec = camera_right * (sensor_width / width);
    const vec3f y_vec = camera_up * (-sensor_height / height);
    const vec3f sensor_pos = sensor_tlc + x_vec * pixel.x() + y_vec * pixel.y();
    const vec3f sensor_dir = sensor_pos - pinhole_pos;

    ray r = { sensor_pos, sensor_dir * (1 / sensor_dir.magnitude()) };

    const double wavelength = chooseRandWavelength();

    spectrum radiance = 0;
    spectrum throughput = 1;

    int dim = 2;

    for (int i = 0; i < max_bounces; ++i)
    {
        const auto scene_hit = world_ptr->nearest_intersection(r);

        const scene_object * const obj = scene_hit.first;
        const float obj_t = scene_hit.second;

        if (!obj)
            break;

        const vec3f surface_pt = r.o + r.d * obj_t;
        const vec3f normal = obj->get_normal(surface_pt);

        radiance += throughput * obj->emission; // Add the emission term from the rendering equation

        const vec3f w_i = -r.d;

        const float n_dot_w_i = dot(normal, w_i);

        if (n_dot_w_i < 0)
            break;

        // HACK
        const float r0 = 0.02f;
        const float fresnel = r0 + (1 - r0) * std::pow(1 - n_dot_w_i, 5.0f);

        const bool use_mirror = rng->nextFloat() < fresnel;

        // Add in estimate of reflected light (attenuate for next bounce's emission):
        vec3f random_dir;
        spectrum weight;
        if (obj->mirror && use_mirror)
        {
            random_dir = normal * (2 * dot(normal, w_i)) - w_i;
            weight = 1;
        }
        else
        {
            // Step 1: Get uniform reflection vector in same hemi sphere as normal

            // From Global Illumination Compendium, eq 33.
            const float r1 = randomized_halton(sample_idx, dim++, pixel_idx); // [0,1)
            const float r2 = randomized_halton(sample_idx, dim++, pixel_idx);
            const float k = 2 * std::sqrt(r2 * (1 - r2));
            const float phi = float_twopi * r1;
            random_dir = vec3f(
                k * std::cos(phi),
                1 - 2 * r2,
                k * std::sin(phi)
            );

            // Flip if in wrong hemisphere.
            const float n_dot_w_o = dot(random_dir, normal);
            if (n_dot_w_o < 0)
                random_dir = -random_dir;
            // Step 2: Weight the throughput by BRDF * cos / p.

            // w = BRDF * cos / p
            // = (albedo/pi) * cos / (1 / 2pi)
            weight = (obj->albedo / float_pi) * std::fabs(n_dot_w_o) * float_twopi;
            throughput *= weight;
        }

        // Time for some Russian roulette!
        const float max_weight = std::max(std::max(weight.x(), weight.y()), weight.z());
        const float rr_live_p = std::min(1.0f, max_weight);
        if (randomized_halton(sample_idx, dim++, pixel_idx) > rr_live_p)
            break;
        else
            throughput *= (1 / rr_live_p);

        // Start tracing from the new point.
        r.o = surface_pt;
        r.d = random_dir;
    }

    return vec3f(radiance.x(), radiance.y(), radiance.z());

    //const auto scene_hit = world_ptr->nearest_intersection(r);
    //const auto obj = scene_hit.first;
    //const auto obj_t = scene_hit.second;
    //if(obj) {
    //    const auto surface_pt = r.o + r.d * obj_t;
    //    const auto norm = obj->get_normal(surface_pt);
    //    const auto intensity = dot(r.d * -1, norm);
    //    return obj->albedo * intensity;
    //}
    //return 0;
}

void ThreadFunc(const int thread_id, pcg32 * rng, const int width, const int height, std::atomic<int>* curr_block, world* const world_container, uint8_t* const pixels)
{
    constexpr int block_size = 16;
    const int num_samples = 1 << 8;

    rng->seed(thread_id + 1337); // Don't pass 0 as rng seed.

    //Initialize blocks
    while(true)
    {
        const int x_blocks = (width  + block_size - 1) / block_size; // Round up the number of blocks.
        const int y_blocks = (height + block_size - 1) / block_size;

        const int block_idx = curr_block->fetch_add(1);

        if (block_idx >= x_blocks * y_blocks)
            return;

        const int block_y  = block_idx / x_blocks;
        const int block_x  = block_idx - x_blocks * block_y; // Equiv to block_idx % x_blocks.
        const int block_x0 = block_x * block_size, block_x1 = std::min(block_x0 + block_size, width);
        const int block_y0 = block_y * block_size, block_y1 = std::min(block_y0 + block_size, height);

        for (int j = block_x0; j < block_x1; ++j)
        for (int i = block_y0; i < block_y1; ++i)
        {
            vec3f sum = 0;
            for (int s = 0; s < num_samples; ++s) {
                const vec3f color_linear = pathtrace(i, j, s, width, height, world_container, rng);
                sum += color_linear;
            }
            if(j < height / 6) {
            const vec3f wavelength_color = spectrum::xyz2rgb(spectrum::wavelength2xyz(wavelength_min + (static_cast<double>(i)/width) * (wavelength_max - wavelength_min)));
            sum += wavelength_color * num_samples;
            }
            const vec3f color_linear = sum / num_samples;
            //const vec3f color_linear = spectrum::xyz2rgb(color_xyz);

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


int main(int argc, char* argv[])
{
	const int width  = 512;
	const int height = 512;
	const int channels = 3;
	const int num_samples = 1 << 8;

    constexpr vec3f origin(0,0,0);

    //std::unique_ptr<sphere>


    const int num_frames = 1;

    for (int frame = 0; frame < num_frames; ++frame)
    {
        std::unique_ptr<sphere> obj1 = std::make_unique<sphere>(
            spectrum{0.85f, 0.25f, 0.05f}, // albedo
            spectrum{0.0f}, // emission
            3.0f, // radius
            origin
            );
        obj1->mirror = true;

        const float a = frame * float_twopi / num_frames;
        std::unique_ptr<sphere> obj2 = std::make_unique<sphere>(
            spectrum(0.85f), // albedo
            spectrum{0.0f}, // emission
            1.0f, // radius
            vec3f{ 0, 5, 0 } + vec3f{ std::cos(a), std::sin(a), 0 }
        );
        obj2->mirror = true;

        const std::unique_ptr<sphere> obj3 = std::make_unique<sphere>(
            spectrum(0), // albedo
            spectrum{1.0f}, // emission
            30.0f, // radius
            vec3f{ 0, 40, 0 }
        );

        std::unique_ptr<sphere> obj4 = std::make_unique<sphere>(
            spectrum(0.85f), // albedo
            spectrum{0.0f}, // emission
            1.0f, // radius
            vec3f{ 0, 4, 0 }
        );
        obj2->mirror = true;

        std::unique_ptr<world> world_container = std::make_unique<world>();
        world_container->objects.push_back(obj1.get());
        world_container->objects.push_back(obj2.get());
        world_container->objects.push_back(obj3.get());

        std::vector<uint8_t> pixels(width * height * channels);

        std::atomic<int> curr_block(0);

        const int num_threads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads;
        std::vector<pcg32> rngs(num_threads);
        threads.reserve(num_threads);

        constexpr int block_size = 16;
        //const auto start_time = std::chrono::
        for(int thread_id = 0; thread_id < num_threads; ++thread_id) {
            // void ThreadFunc(const int thread_id, pcg32 * rng, const int width, const int height, std::atomic<int>* curr_block, world* const world_container, uint8_t* const pixels)

            threads.push_back(std::thread(ThreadFunc, thread_id, &rngs[thread_id], width, height, &curr_block, world_container.get(), &pixels[0]));
        }

        for(auto& t : threads){
            t.join();
        }


        char filename[256];
        sprintf(filename, "frame%.4d.png", frame);
        stbi_write_png(filename, width, height, channels, &pixels[0], width * channels);
        printf("wrote %s\n", filename);
    }

	return 0;
}
