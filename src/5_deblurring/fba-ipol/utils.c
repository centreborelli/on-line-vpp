// Copyright (C) 2016, Jérémy Anger <jeremy.anger@cmla.ens-cachan.fr>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "utils.h"

#define PIXEL3(img, x, y, z) ((img).data[(z) + (img).d * ((x) + (img).w * (y))])
#define PIXEL2(img, x, y) ((img).data[(x) + (img).w * (y)])
#define PIXEL1(img, i) ((img).data[(i)])
// trick from http://stackoverflow.com/a/11763277
#define GET_MACRO(_1,_2,_3,_4,NAME,...) NAME
#define PIXEL(...) GET_MACRO(__VA_ARGS__, PIXEL3, PIXEL2, PIXEL1, NOP)(__VA_ARGS__)

image_float_t new_image_float(int w, int h, int d)
{
    return (image_float_t) {
        .w=w, .h=h, .d=d, .size=w*h*d*sizeof(float), .type="float",
            .data=(float*) malloc(sizeof(float)*w*h*d)
    };
}
image_double_t new_image_double(int w, int h, int d)
{
    return (image_double_t) {
        .w=w, .h=h, .d=d, .size=w*h*d*sizeof(double), .type="double",
            .data=(double*) malloc(sizeof(double)*w*h*d)
    };
}
image_complex_t new_image_complex(int w, int h, int d)
{
    return (image_complex_t) {
        .w=w, .h=h, .d=d, .size=w*h*d*sizeof(double complex), .type="complex",
            .data=(double complex*) fftw_malloc(sizeof(double complex)*w*h*d)
    };
}

image_float_t new_image_float_data(int w, int h, int d, float *data)
{
    return (image_float_t) {
        .w=w, .h=h, .d=d, .size=w*h*d*sizeof(float), .type="float",
        .data=data,
    };
}

void set_value(image_double_t img, double value)
{
    for (int j = 0; j < img.w*img.h*img.d; j++)
        PIXEL(img, j) = value;
}

void set_valuef(image_float_t img, float value)
{
    for (int j = 0; j < img.w*img.h*img.d; j++)
        PIXEL(img, j) = value;
}

void set_valuec(image_complex_t img, double complex value)
{
    for (int j = 0; j < img.w*img.h*img.d; j++)
        PIXEL(img, j) = value;
}

void copyf(image_float_t out, image_float_t in)
{
    assert(out.w == in.w);
    assert(out.h == in.h);
    assert(out.d == in.d);

    for (int j = 0; j < out.w*out.h*out.d; j++)
        PIXEL(out, j) = PIXEL(in, j);
}

void convert_f2d(image_double_t out, image_float_t in)
{
    assert(out.w == in.w);
    assert(out.h == in.h);
    assert(out.d == in.d);

    for (int j = 0; j < out.w*out.h*out.d; j++)
        PIXEL(out, j) = PIXEL(in, j);
}

void convert_d2f(image_float_t out, image_double_t in)
{
    assert(out.w == in.w);
    assert(out.h == in.h);
    assert(out.d == in.d);

    for (int j = 0; j < out.w*out.h*out.d; j++)
        PIXEL(out, j) = PIXEL(in, j);
}

void padding(image_double_t out, image_float_t in)
{
    assert(out.w >= in.w);
    assert(out.h >= in.h);
    assert(out.d == in.d);

    set_value(out, 0.);

    for (int y = 0; y < out.h; y++)
    for (int x = 0; x < out.w; x++)
    {
        int xx = x;
        int yy = y;
        if (yy >= in.h)
            yy = in.h - (yy - in.h + 1);
        if (xx >= in.w)
            xx = in.w - (xx - in.w + 1);
        for (int dd = 0; dd < in.d; dd++)
            PIXEL(out, x, y, dd) = PIXEL(in, xx, yy, dd);
    }
}

void paddingf(image_float_t out, image_float_t in)
{
    assert(out.w >= in.w);
    assert(out.h >= in.h);
    assert(out.d == in.d);

    set_valuef(out, 0.f);

    for (int y = 0; y < out.h; y++)
    for (int x = 0; x < out.w; x++)
    {
        int xx = x;
        int yy = y;
        if (yy >= in.h)
            yy = in.h - (yy - in.h + 1);
        if (xx >= in.w)
            xx = in.w - (xx - in.w + 1);
        for (int dd = 0; dd < in.d; dd++)
            PIXEL(out, x, y, dd) = PIXEL(in, xx, yy, dd);
    }
}

void extract(image_double_t patch, image_double_t img, int x, int y)
{
    assert(patch.d == img.d);
    assert(patch.w + x <= img.w);
    assert(patch.h + y <= img.h);

    int i = 0;
    for (int yy = y; yy < patch.h + y; yy++)
    for (int xx = x; xx < patch.w + x; xx++)
    for (int dd = 0; dd < img.d; dd++)
        PIXEL(patch, i++) = PIXEL(img, xx, yy, dd);
}

void fft(fftw_plan plan, image_complex_t ft, image_double_t img, image_complex_t buffer)
{
    assert(ft.d == 3);
    assert(img.d == 3);

    for (int d = 0; d < 3; d++)
    {
        for (int j = 0; j < img.w * img.h; j++)
            PIXEL(buffer, j) = PIXEL(img, j*3+d); // real to complex
        fftw_execute_dft(plan, buffer.data, buffer.data);
        for (int j = 0; j < img.w * img.h; j++)
            PIXEL(ft, j*3+d) = PIXEL(buffer, j);
    }
}

void ifft(fftw_plan plan, image_double_t img, image_complex_t ft, image_complex_t buffer)
{
    assert(ft.d == 3);
    assert(img.d == 3);

    float norm = ft.w * ft.h;
    for (int d = 0; d < 3; d++)
    {
        // normalize matlab style
        for (int j = 0; j < img.w * img.h; j++)
            PIXEL(buffer, j) = PIXEL(ft, j*3+d) / norm;
        fftw_execute_dft(plan, buffer.data, buffer.data);
        for (int j = 0; j < img.w * img.h; j++)
            PIXEL(img, j*3+d) = creal(PIXEL(buffer, j));
    }
}

void fftshift(image_double_t img)
{
    assert(img.d == 1);
    assert(img.w % 2 == 0);
    assert(img.h % 2 == 0);

    int halfw = img.w / 2;
    int halfh = img.h / 2;

    for (int y = 0; y < halfh; y++)
    for (int x = 0; x < halfw; x++)
    {
        double tmp;

        // swap img(x, y) with img(x + w/2, y + h/2)
        tmp = PIXEL(img, x, y);
        PIXEL(img, x, y) = PIXEL(img, x + halfw, y + halfh);
        PIXEL(img, x + halfw, y + halfh) = tmp;

        // swap img(x + w/2, y) with img(x, y + h/2)
        tmp = PIXEL(img, x + halfw, y);
        PIXEL(img, x + halfw, y) = PIXEL(img, x, y + halfh);
        PIXEL(img, x, y + halfh) = tmp;
    }
}

void accumulate(image_double_t out, image_double_t in, int x, int y)
{
    assert(out.d == in.d);

    int i = 0;
    for (int yy = y; yy < in.h + y; yy++)
    for (int xx = x; xx < in.w + x; xx++)
    for (int dd = 0; dd < out.d; dd++)
        PIXEL(out, xx, yy, dd) += PIXEL(in, i++);
}

void crop(image_float_t out, image_double_t in)
{
    assert(out.w <= in.w);
    assert(out.h <= in.h);
    assert(out.d == in.d);

    for (int y = 0; y < out.h; y++)
    for (int x = 0; x < out.w; x++)
    for (int dd = 0; dd < out.d; dd++)
        PIXEL(out, x, y, dd) = PIXEL(in, x, y, dd);
}

void cropf(image_float_t out, image_float_t in)
{
    assert(out.w <= in.w);
    assert(out.h <= in.h);
    assert(out.d == in.d);

    for (int y = 0; y < out.h; y++)
    for (int x = 0; x < out.w; x++)
    for (int dd = 0; dd < out.d; dd++)
        PIXEL(out, x, y, dd) = PIXEL(in, x, y, dd);
}

void rgb2greyf(image_float_t grey, image_float_t rgb)
{
    assert(rgb.d == 3);
    assert(grey.d == 1);

    for (int i = 0; i < grey.w * grey.h; i++)
        PIXEL(grey, i) = (PIXEL(rgb, i*3) + PIXEL(rgb, i*3+1) + PIXEL(rgb, i*3+2)) / 3.f;
}

void linear_combination(image_float_t out, image_float_t in1, image_float_t in2, image_float_t mask)
{
    assert(in1.w == out.w);
    assert(in1.h == out.h);
    assert(in1.d == out.d);
    assert(in2.w == out.w);
    assert(in2.h == out.h);
    assert(in2.d == out.d);
    assert(mask.w == out.w);
    assert(mask.h == out.h);
    assert(mask.d == 1);

    for (int i = 0; i < out.w*out.h; i++)
    {
        float m = PIXEL(mask, i);
        for (int dd = 0; dd < out.d; dd++)
        {
            int index = i*out.d+dd;
            PIXEL(out, index) = m * PIXEL(in1, index) + (1 - m) * PIXEL(in2, index);
        }
    }
}

#include "bicubic.c"
void upsamplef(image_float_t out, image_float_t in)
{
    assert(out.d == in.d);
    assert(out.w >= in.w);
    assert(out.h >= in.h);

    float rx = (float) in.w / out.w;
    float ry = (float) in.h / out.h;

    for (int y = 0; y < out.h; y++)
    for (int x = 0; x < out.w; x++)
        bicubic_interpolation(&PIXEL(out, x, y, 0), in.data, in.w, in.h, in.d, x * rx, y * ry);
}

void downsamplef(image_float_t out, image_float_t in)
{
    assert(out.d == in.d);
    assert(out.w <= in.w);
    assert(out.h <= in.h);
    assert(in.w % out.w == 0);
    assert(in.h % out.h == 0);

    int rx = in.w / out.w;
    int ry = in.h / out.h;

    for (int y = 0; y < out.h; y++)
    for (int x = 0; x < out.w; x++)
    for (int dd = 0; dd < out.d; dd++)
    {
        double sum = 0.;
        int num = 0;
        for (int yy = 0; yy < ry; yy++)
        for (int xx = 0; xx < rx; xx++)
            if (x*rx+xx < in.w && y*ry+yy < in.h)
            {
                sum += PIXEL(in, xx+rx*x, yy+ry*y, dd);
                num++;
            }
        PIXEL(out, x, y, dd) = sum / num;
    }
}

#include "erosion.c"
void erodef_disk5(image_float_t out, image_float_t in)
{
    assert(in.d == 1);
    assert(out.d == 1);
    assert(in.w == out.w);
    assert(in.h == out.h);

    static int* disk = 0;
    if (!disk)
        disk = build_disk(5);
    morsi_erosion(out.data, in.data, in.w, in.h, disk);
}

