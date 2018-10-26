// Copyright (C) 2016, Jérémy Anger <jeremy.anger@cmla.ens-cachan.fr>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include <complex.h>
#include <fftw3.h>

#include "gaussian_conv_vyv.h"
#include "utils.h"
#include "fba.h"

struct fba_process_inputs {
    int x;
    int y;
    int n;
    float p;
    image_double_t* buffer_padded;
    image_double_t window;
};

// to store intermediate results
struct fba_process_buffers {
    image_double_t average_magnitude;
    image_double_t weights_accum;

    image_complex_t ft_accum;
    image_complex_t ft_patch;
    image_complex_t ft_buffer;

    vyv_coeffs blur_coeffs;
    fftw_plan fft_forward_plan;
    fftw_plan fft_backward_plan;
};

void compute_patch_from_temporal_window(image_double_t patch, struct fba_process_inputs* pi, struct fba_process_buffers* pb)
{
    int W = patch.w;
    int d = patch.d;

    set_valuec(pb->ft_accum, 0.);
    set_value(pb->weights_accum, 0.);

    for (int m = 0; m < pi->n; m++)
    {
        // extract the patch from the image at index 'm'
        extract(patch, pi->buffer_padded[m], pi->x, pi->y);

        // compute the fft of the patch
        fft(pb->fft_forward_plan, pb->ft_patch, patch, pb->ft_buffer);

        // compute the average magnitude of the Fourier transform
        for (int i = 0; i < W * W; i++)
        {
            double magnitude = 0.;
            for (int dd = 0; dd < d; dd++)
                magnitude += cabs(pb->ft_patch.data[i*d+dd]);
            pb->average_magnitude.data[i] = magnitude / d;
        }

        // blur the average magnitude
        // NOTE: fftshift is used to group up low frequencies together
        //       because vyv_gaussian_conv_image is not a periodic convolution
        fftshift(pb->average_magnitude);
        vyv_gaussian_conv_image(pb->blur_coeffs, pb->average_magnitude.data, pb->average_magnitude.data, W, W, 1);
        fftshift(pb->average_magnitude);

        // accumulate the weighted Fourier transform
        for (int j = 0; j < W*W; j++)
        {
            double weight = pow(pb->average_magnitude.data[j], pi->p);
            for (int dd = 0; dd < d; dd++)
                pb->ft_accum.data[j*d+dd] += pb->ft_patch.data[j*d+dd] * weight;
            pb->weights_accum.data[j] += weight;
        }
    }

    // divide by the weights sum
    for (int j = 0; j < W*W; j++)
    for (int dd = 0; dd < d; dd++)
        pb->ft_accum.data[j*d+dd] /= pb->weights_accum.data[j] + 1e-8;

    // compute the inverse fft
    ifft(pb->fft_backward_plan, patch, pb->ft_accum, pb->ft_buffer);

    // weight the result using the window
    for (int j = 0; j < W*W; j++)
    {
        double v = pi->window.data[j];
        for (int dd = 0; dd < d; dd++)
            patch.data[j*d + dd] *= v;
    }
}

void fba(image_float_t out, image_float_t* inputs, int W, float p, int n)
{
    assert(inputs[0].d == 3);

    const int w = out.w;
    const int h = out.h;
    const int d = out.d;
    // pad the input to have an integer number of tiles
    const int pw = W/2 * (w / (W/2) + 1);
    const int ph = W/2 * (h / (W/2) + 1);

    // precompute the blurring coefficients (to smooth the weights)
    double sigma = W / 50.0;
    vyv_coeffs blur_coeffs;
    vyv_precomp(&blur_coeffs, sigma, 3, 0.01);

    // prepare the fftw plans
    fftw_plan fft_forward_plan, fft_backward_plan;
#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
    {
        image_complex_t buffer = new_image_complex(W, W, 1);
        fft_forward_plan = fftw_plan_dft_2d(W, W, buffer.data, buffer.data, FFTW_FORWARD, FFTW_MEASURE);
        fft_backward_plan = fftw_plan_dft_2d(W, W, buffer.data, buffer.data, FFTW_BACKWARD, FFTW_MEASURE);
        fftw_free(buffer.data);
    }

    image_double_t buffer_padded[n];
    image_double_t output_padded = new_image_double(pw, ph, d);
    image_double_t windowing_accum = new_image_double(pw, ph, 1);
    // set the window to a constant image
    image_double_t window = new_image_double(W, W, 1);
    set_value(window, 1.);

    // pad the input
    for (int m = 0; m < n; m++)
    {
        buffer_padded[m] = new_image_double(pw, ph, d);
        padding(buffer_padded[m], inputs[m]);
    }

    // set the accumulators to 0
    set_value(output_padded, 0.);
    set_value(windowing_accum, 0.);

    image_double_t patch = new_image_double(W, W, d);

    // allocate temporary images
    struct fba_process_buffers process_buffers;
    process_buffers.ft_patch = new_image_complex(W, W, 3);
    process_buffers.ft_accum = new_image_complex(W, W, d);
    process_buffers.average_magnitude = new_image_double(W, W, 1);
    process_buffers.weights_accum = new_image_double(W, W, 1);
    process_buffers.ft_buffer = new_image_complex(W, W, 1);

    process_buffers.blur_coeffs = blur_coeffs;
    process_buffers.fft_forward_plan = fft_forward_plan;
    process_buffers.fft_backward_plan = fft_backward_plan;

    // process each tile and accumulate the result into output_padded
    for (int y = 0; y < h - W/2; y += W / 2)
    {
        for (int x = 0; x < w - W/2; x += W / 2)
        {
            struct fba_process_inputs process_inputs;
            process_inputs.x = x;
            process_inputs.y = y;
            process_inputs.n = n;
            process_inputs.p = p;
            process_inputs.buffer_padded = buffer_padded;
            process_inputs.window = window;
            compute_patch_from_temporal_window(patch, &process_inputs, &process_buffers);

            // accumulate the result and the window
            accumulate(output_padded, patch, x, y);
            accumulate(windowing_accum, window, x, y);
        }
    }

    // divide the result by the windowing weights
    for (int j = 0; j < pw*ph; j++)
    {
        double v = windowing_accum.data[j];
        for (int dd = 0; dd < d; dd++)
            output_padded.data[j*d + dd] /= v;
    }

    // extract the interesting part of the image
    crop(out, output_padded);

    free(patch.data);
    free(process_buffers.ft_patch.data);
    free(process_buffers.ft_accum.data);
    free(process_buffers.average_magnitude.data);
    free(process_buffers.weights_accum.data);
    fftw_free(process_buffers.ft_buffer.data);

    free(window.data);
    free(output_padded.data);
    free(windowing_accum.data);

    for (int m = 0; m < n; m++)
        free(buffer_padded[m].data);

#ifdef _OPENMP
#pragma omp critical (fftw)
#endif
    {
        fftw_destroy_plan(fft_forward_plan);
        fftw_destroy_plan(fft_backward_plan);
    }
}

