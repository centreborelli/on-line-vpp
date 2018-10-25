// Copyright (C) 2016, Jérémy Anger <jeremy.anger@cmla.ens-cachan.fr>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "gaussian_conv_vyv.h"
#include "fba.h"
#include "consistent_registration.h"

void Dual_TVL1_optic_flow_multiscale(
		float *I0,           // source image
		float *I1,           // target image
		float *u1,           // x component of the optical flow
		float *u2,           // y component of the optical flow
		const int   nxx,     // image width
		const int   nyy,     // image height
		const float tau,     // time step
		const float lambda,  // weight parameter for the data term
		const float theta,   // weight parameter for (u - v)²
		const int   nscales, // number of scales
		const float zfactor, // factor for building the image piramid
		const int   warps,   // number of warpings per scale
		const float epsilon, // tolerance for numerical convergence
		const bool  verbose,  // enable/disable the verbose mode
        const int max_iterations
);

void bicubic_interpolation_warp(
	const float *input,     // image to be warped
	const float *u,         // x component of the vector field
	const float *v,         // y component of the vector field
	float       *output,    // image warped with bicubic interpolation
	const int    nx,        // image width
	const int    ny,        // image height
	bool         border_out // if true, put zeros outside the region
);

float bicubic_interpolation_at(
	const float *input, //image to be interpolated
	const float  uu,    //x component of the vector field
	const float  vv,    //y component of the vector field
	const int    nx,    //image width
	const int    ny,    //image height
	bool         border_out //if true, return zero outside the region
);

void consistent_registration(image_float_t* outputs, image_float_t* inputs, int n, int ref, float eps, int downsampling)
{
    const int w = inputs[0].w;
    const int h = inputs[0].h;
    const int d = inputs[0].d;
    assert(d == 3);

    // compute the padded size and downscaled size
    const int ww = w + (downsampling - w % downsampling);
    const int hh = h + (downsampling - h % downsampling);
    const int sw = ww / downsampling;
    const int sh = hh / downsampling;

    image_float_t buffer[n];
    image_float_t buffer_grey[n];
    image_float_t buffer_small[n];

    for (int m = 0; m < n; m++)
    {
        buffer[m] = new_image_float(ww, hh, d);
        buffer_grey[m] = new_image_float(ww, hh, 1);
        buffer_small[m] = new_image_float(sw, sh, 1);
    }

    // pad, convert to gray and downsample the inputs
    for (int m = 0; m < n; m++)
    {
        paddingf(buffer[m], inputs[m]);
        rgb2greyf(buffer_grey[m], buffer[m]);
        downsamplef(buffer_small[m], buffer_grey[m]);
    }

    for (int m = 0; m < n; m++)
        set_valuef(outputs[m], 0.f);

    // compute blurring coefficients (for the consistency map)
    float sigma = 5.f;
    vyv_coeffs blur_coeffs;
    vyv_precomp(&blur_coeffs, sigma, 3, 0.01);

    // allocate auxiliary images
    image_float_t fflowx = new_image_float(ww, hh, 1);
    image_float_t fflowy = new_image_float(ww, hh, 1);
    image_float_t bflowx = new_image_float(ww, hh, 1);
    image_float_t bflowy = new_image_float(ww, hh, 1);
    image_float_t fflowx_small = new_image_float(sw, sh, 1);
    image_float_t fflowy_small = new_image_float(sw, sh, 1);
    image_float_t bflowx_small = new_image_float(sw, sh, 1);
    image_float_t bflowy_small = new_image_float(sw, sh, 1);
    image_float_t cmap = new_image_float(ww, hh, 1);
    image_float_t registered = new_image_float(ww, hh, d);
    image_float_t buf_grey = new_image_float(ww, hh, 1);
    image_float_t registered_grey = new_image_float(ww, hh, 1);
    image_float_t cmap_eroded = new_image_float(ww, hh, 1);
    image_double_t cmap_eroded_double = new_image_double(ww, hh, 1);
    image_double_t cmap_blurred_double = new_image_double(ww, hh, 1);

    for (int m = 0; m < n; m++)
    {
        // if the current image isn't the reference, we register it
        if (m != ref)
        {
            // forward and backward optical flow estimation (on downscaled images)
            Dual_TVL1_optic_flow_multiscale(buffer_small[ref].data, buffer_small[m].data, fflowx_small.data, fflowy_small.data,
                                            sw, sh, 0.25, 0.15, 0.3, 5, 0.5, 3, 0.001, 0, 10);
            Dual_TVL1_optic_flow_multiscale(buffer_small[m].data, buffer_small[ref].data, bflowx_small.data, bflowy_small.data,
                                            sw, sh, 0.25, 0.15, 0.3, 5, 0.5, 3, 0.001, 0, 10);

            // upscale the flow (bicubic interpolation) and multiply by downsampling to keep the correct scale
            upsamplef(fflowx, fflowx_small);
            upsamplef(fflowy, fflowy_small);
            upsamplef(bflowx, bflowx_small);
            upsamplef(bflowy, bflowy_small);
            for (int i = 0; i < ww*hh; i++)
            {
                fflowx.data[i] *= (float) downsampling;
                fflowy.data[i] *= (float) downsampling;
                bflowx.data[i] *= (float) downsampling;
                bflowy.data[i] *= (float) downsampling;
            }

            // compute the registration error
            for (int y = 0; y < hh; y++)
            for (int x = 0; x < ww; x++)
            {
                // xx,yy = displaced x,y by forward flow
                float xx = x + bicubic_interpolation_at(fflowx.data, x, y, ww, hh, 0);
                float yy = y + bicubic_interpolation_at(fflowy.data, x, y, ww, hh, 0);
                // xxx,yyy = displaced xx,yy by backflow
                float xxx = xx + bicubic_interpolation_at(bflowx.data, xx, yy, ww, hh, 0);
                float yyy = yy + bicubic_interpolation_at(bflowy.data, xx, yy, ww, hh, 0);
                // if xxx,yyy is different than x,y, the pixel is marked as inconsistent
                cmap.data[x+y*ww] = hypotf(xxx - x, yyy - y) < eps ? 1.f : 0.f;
            }

            // erode and blur the consistency map
            erodef_disk5(cmap_eroded, cmap);
            convert_f2d(cmap_eroded_double, cmap_eroded);
            vyv_gaussian_conv_image(blur_coeffs, cmap_blurred_double.data, cmap_eroded_double.data, ww, hh, 1 /* d */);
            convert_d2f(cmap, cmap_blurred_double);

            // warp the frame, channel by channel
            for (int dd = 0; dd < d; dd++)
            {
                for (int i = 0; i < ww*hh; i++)
                    buf_grey.data[i] = buffer[m].data[i*d+dd];
                // interpolate using the forward flow
                bicubic_interpolation_warp(buf_grey.data, fflowx.data, fflowy.data, registered_grey.data, ww, hh, true);
                for (int i = 0; i < ww*hh; i++)
                {
                    registered.data[i*d+dd] = registered_grey.data[i];
                    // if the flow led outside the image, indicate that the area is inconsistency
                    if (registered.data[i*d+dd] == 0.f)
                        cmap.data[i] = 0.f;
                }
            }

            // replace warped frame by buffer[registration] if the consistency is low
            linear_combination(registered, registered, buffer[ref], cmap);

            // crop the result
            cropf(outputs[m], registered);
        }
        else
        {
            // if the image is already the reference, just copy it instead of registrating it
            copyf(outputs[m], inputs[m]);
        }
    }

    free(fflowx.data);
    free(fflowy.data);
    free(bflowx.data);
    free(bflowy.data);
    free(fflowx_small.data);
    free(fflowy_small.data);
    free(bflowx_small.data);
    free(bflowy_small.data);
    free(cmap.data);
    free(registered.data);
    free(buf_grey.data);
    free(registered_grey.data);
    free(cmap_eroded.data);
    free(cmap_eroded_double.data);
    free(cmap_blurred_double.data);

    for (int m = 0; m < n; m++)
    {
        free(buffer[m].data);
        free(buffer_grey[m].data);
        free(buffer_small[m].data);
    }
}

