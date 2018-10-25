#ifndef UTILS
#define UTILS

#include <complex.h>
#include <fftw3.h>

typedef struct {
    int w, h, d;
    size_t size;
    const char* type;
    float* data;
} image_float_t;

typedef struct {
    int w, h, d;
    size_t size;
    const char* type;
    double* data;
} image_double_t;

typedef struct {
    int w, h, d;
    size_t size;
    const char* type;
    double complex* data;
} image_complex_t;

// allocate a new image
image_float_t new_image_float(int w, int h, int d);
image_double_t new_image_double(int w, int h, int d);
image_complex_t new_image_complex(int w, int h, int d);

// create an image from existing data
image_float_t new_image_float_data(int w, int h, int d, float *data);

// set an image to a constant value
void set_value(image_double_t img, double value);
void set_valuef(image_float_t img, float value);
void set_valuec(image_complex_t img, double complex value);

// copy an image to another
void copyf(image_float_t out, image_float_t in);

// convert an image from float to double
void convert_f2d(image_double_t out, image_float_t in);
// convert an image from double to float
void convert_d2f(image_float_t out, image_double_t in);

// pad an image using symmetric boundary
void padding(image_double_t out, image_float_t in);
void paddingf(image_float_t out, image_float_t in);

// extract a patch from an image
void extract(image_double_t patch, image_double_t img, int x, int y);

// compute the fast fourier transform and its inverse
void fft(fftw_plan plan, image_complex_t ft, image_double_t img, image_complex_t buffer);
void ifft(fftw_plan plan, image_double_t img, image_complex_t ft, image_complex_t buffer);

// circular shift of an halfsize
void fftshift(image_double_t img);

// accumulate a patch at a given position
void accumulate(image_double_t out, image_double_t in, int x, int y);

// crop an image
void crop(image_float_t out, image_double_t in);
void cropf(image_float_t out, image_float_t in);

// convert an RGB image to greyscale
void rgb2greyf(image_float_t grey, image_float_t rgb);

// linear combination of two image (out = mask*in1 + (1-mask)*in2)
void linear_combination(image_float_t out, image_float_t in1, image_float_t in2, image_float_t mask);

// upsample or downsample an image
void upsamplef(image_float_t out, image_float_t in);
void downsamplef(image_float_t out, image_float_t in);

// erode an image by a disk of radius 5
void erodef_disk5(image_float_t out, image_float_t in);

#endif

