/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file utilities.cpp
 * @brief Utilities functions.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "Utilities.h"

#include <math.h>
#include <omp.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>

#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3

using namespace std;

/**
 * @brief Convenient function to use the sort function provided by the vector library.
 **/
bool comparaisonFirst(
	const pair<float, unsigned> &i_pair1
,	const pair<float, unsigned> &i_pair2
){
	return i_pair1.first < i_pair2.first;
}

/**
 * @brief Check if a number is a power of 2
 **/
bool power_of_2(
    const unsigned n
){
    if (n == 0)
        return false;

    if (n == 1)
        return true;

    if (n % 2 == 0)
        return power_of_2(n / 2);
    else
        return false;
}

/**
 * @brief Look for the closest power of 2 number
 *
 * @param n: number
 *
 * @return the closest power of 2 lower or equal to n
 **/
int closest_power_of_2(
    const unsigned n
){
    unsigned r = 1;
    while (r * 2 <= n)
        r *= 2;

    return r;
}

/**
 * @brief Initialize a 2D fftwf_plan with some parameters
 *
 * @param plan: fftwf_plan to allocate;
 * @param N: size of the patch to apply the 2D transform;
 * @param kind: forward or backward;
 * @param nb: number of 2D transform which will be processed.
 *
 * @return none.
 **/
void allocate_plan_2d(
    fftwf_plan* plan
,   const unsigned N
,   const fftwf_r2r_kind kind
,   const unsigned nb
){
    int            nb_table[2]   = {N, N};
    int            nembed[2]     = {N, N};
    fftwf_r2r_kind kind_table[2] = {kind, kind};

    float* vec = (float*) fftwf_malloc(N * N * nb * sizeof(float));
    (*plan) = fftwf_plan_many_r2r(2, nb_table, nb, vec, nembed, 1, N * N, vec,
                                  nembed, 1, N * N, kind_table, FFTW_ESTIMATE);

    fftwf_free(vec);
}

/**
 * @brief Initialize a 1D fftwf_plan with some parameters
 *
 * @param plan: fftwf_plan to allocate;
 * @param N: size of the vector to apply the 1D transform;
 * @param kind: forward or backward;
 * @param nb: number of 1D transform which will be processed.
 *
 * @return none.
 **/
void allocate_plan_1d(
    fftwf_plan* plan
,   const unsigned N
,   const fftwf_r2r_kind kind
,   const unsigned nb
){
    int nb_table[1] = {N};
    int nembed[1]   = {N * nb};
    fftwf_r2r_kind kind_table[1] = {kind};

    float* vec = (float*) fftwf_malloc(N * nb * sizeof(float));
    (*plan) = fftwf_plan_many_r2r(1, nb_table, nb, vec, nembed, 1, N, vec,
                                  nembed, 1, N, kind_table, FFTW_ESTIMATE);
    fftwf_free(vec);
}

/**
 * @brief Initialize a set of indices.
 *
 * @param ind_set: will contain the set of indices;
 * @param max_size: indices can't go over this size;
 * @param N : boundary;
 * @param step: step between two indices.
 *
 * @return none.
 **/
void ind_initialize(
    vector<unsigned> &ind_set
,   const unsigned beginning 
,   const unsigned end
,   const unsigned step
){
    ind_set.clear();
    unsigned ind = beginning;
    while (ind <= end)
    {
        ind_set.push_back(ind);
        ind += step;
    }
    if (ind_set.back() < end)
        ind_set.push_back(end);
}

void ind_initialize2(
    vector<unsigned> &ind_set
,   const unsigned max_size
,   const unsigned N
,   const unsigned step
){
    ind_set.clear();
    unsigned ind = N;
    while (ind < max_size - N)
    {
        ind_set.push_back(ind);
        ind += step;
    }
    if (ind_set.back() < max_size - N - 1)
        ind_set.push_back(max_size - N - 1);
}


/**
 * @brief For convenience. Estimate the size of the ind_set vector built
 *        with the function ind_initialize().
 *
 * @return size of ind_set vector built in ind_initialize().
 **/
unsigned ind_size(
    const unsigned beginning 
,   const unsigned end
,   const unsigned step
){
    unsigned ind = beginning;
    unsigned k = 0;
    while (ind <= end)
    {
        k++;
        ind += step;
    }
    if (ind - step < end)
        k++;

    return k;
}

/**
 * Adapted from the code from Marc Lebrun
 * @brief Transform the color space of an image, from RGB to OPP, or vice-versa.
 *
 * @param io_im: image on which the transform will be applied;
 * @param width, height, chnls: size of io_im;
 * @param p_isForward: if true, go from RGB to OPP, otherwise go from OPP to RGB.
 *
 * @return none.
 **/
void transformColorSpace(
	float* io_im
,   const int width
,   const int height
,   const int chnls
,	const bool p_isForward
){
	//! If the image as only one channel, do nothing
	if (chnls == 1) return;

	//! Initialization
	const unsigned wh     = width * height;
	float* imTmp = (float*) malloc(wh * chnls * sizeof * imTmp);

	//! RGB to OPP
	if (p_isForward) {
		if (chnls == 3) {
			const unsigned red   = 0;
			const unsigned green = wh;
			const unsigned blue  = wh * 2;
			const float a = 1.f / sqrtf(3.f);
			const float b = 1.f / sqrtf(2.f);
			const float c = 2.f * a * sqrtf(2.f);

			for (unsigned k = 0; k < wh; k++) {
				imTmp[k + red  ] = a * (io_im[k + red] + io_im[k + green] + io_im[k + blue]);
				imTmp[k + green] = b * (io_im[k + red] - io_im[k + blue]);
				imTmp[k + blue ] = c * (0.25f * io_im[k + red ] - 0.5f * io_im[k + green]
				                      + 0.25f * io_im[k + blue]);
			}
		}
		else { //! chnls == 4
			const unsigned Gr = 0;
			const unsigned R  = wh;
			const unsigned B  = wh * 2;
			const unsigned Gb = wh * 3;
			const float a = 0.5f;
			const float b = 1.f / sqrtf(2.f);

			for (unsigned k = 0; k < wh; k++) {
				imTmp[k + Gr] = a * ( io_im[k + Gr] + io_im[k + R ] +
				                      io_im[k + B ] + io_im[k + Gb]);
				imTmp[k + R ] = b * ( io_im[k + R ] - io_im[k + B ]);
				imTmp[k + B ] = a * (-io_im[k + Gr] + io_im[k + R ] +
				                      io_im[k + B ] - io_im[k + Gb]);
				imTmp[k + Gb] = b * (-io_im[k + Gr] + io_im[k + Gb]);
			}
		}
	}
	//! OPP to RGB
	else {
		if (chnls == 3) {
			const unsigned red   = 0;
			const unsigned green = wh;
			const unsigned blue  = wh * 2;
			const float a = 1.f / sqrtf(3.f);
			const float b = 1.f / sqrtf(2.f);
			const float c = a / b;

			for (unsigned k = 0; k < wh; k++) {
				//! R channel
				imTmp[k + red  ] = a * io_im[k + red] + b * io_im[k + green]
				                               + c * 0.5f * io_im[k + blue];
				//! G channel
				imTmp[k + green] = a * io_im[k + red] - c * io_im[k + blue];

				//! B channel
				imTmp[k + blue ] = a * io_im[k + red] - b * io_im[k + green]
				                               + c * 0.5f * io_im[k + blue];
			}
		}
		else {	//! chnls == 4
			const unsigned Gr = 0;
			const unsigned R  = wh;
			const unsigned B  = wh * 2;
			const unsigned Gb = wh * 3;
			const float a = 0.5f;
			const float b = 1.f / sqrtf(2.f);
			for (unsigned k = 0; k < wh; k++) {
				imTmp[k + Gr] = a * io_im[k + Gr] - a * io_im[k + B] - b * io_im[k + Gb];
				imTmp[k + R ] = a * io_im[k + Gr] + b * io_im[k + R] + a * io_im[k + B];
				imTmp[k + B ] = a * io_im[k + Gr] - b * io_im[k + R] + a * io_im[k + B];
				imTmp[k + Gb] = a * io_im[k + Gr] - a * io_im[k + B] + b * io_im[k + Gb];
			}
		}
	}

    for(int i = 0; i < wh*chnls; ++i)
        io_im[i] = imTmp[i];

    free(imTmp);
}
