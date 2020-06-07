//
// Created by Craig Scott on 6/4/20.
//

#include "vector.h"
#include "stb_write_impl.h"

using byte = unsigned char;
int main(int argc, char *argv[])
{
    int width, height, channels;
    width = 200;
    height = 200;

    channels = 3;
    uint8_t* pixels = new uint8_t[width*height*channels];

    int index = 0;
    for (int j = 0; j < height; ++j)
    {
        for (int i = 0; i < width; ++i)
        {

            float circle = vec2f(width/2 - i,height/2 - j).magnitude();
            float r = circle > 30.0f ? 0.0f : 1.0f ;
            float g = r;
            float b = r;
            int ir = int(255.99 * r);
            int ig = int(255.99 * g);
            int ib = int(255.99 * b);

            pixels[index++] = ir;
            pixels[index++] = ig;
            pixels[index++] = ib;
        }
    }

    stbi_write_png("stbpng.png", width, height, channels, pixels, width * channels);
    delete [] pixels;

    return 0;
}
