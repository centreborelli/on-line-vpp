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
#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include <vector>
#include <fftw3.h>

/**
 * @brief Convenient function to use the sort function provided by the vector library.
 **/
bool comparaisonFirst(
	const std::pair<float, unsigned> &i_pair1
,	const std::pair<float, unsigned> &i_pair2
);

//! Check if a number is a power of 2
bool power_of_2(
    const unsigned n
);

//! Look for the closest power of 2 number
int closest_power_of_2(
    const unsigned n
);

//! Initialize a 2D fftwf_plan with some parameters
void allocate_plan_2d(
    fftwf_plan* plan
,   const unsigned N
,   const fftwf_r2r_kind kind
,   const unsigned nb
);

//! Initialize a 1D fftwf_plan with some parameters
void allocate_plan_1d(
    fftwf_plan* plan
,   const unsigned N
,   const fftwf_r2r_kind kind
,   const unsigned nb
);

//! Initialize a set of indices
void ind_initialize(
    std::vector<unsigned> &ind_set
,   const unsigned beginning
,   const unsigned end
,   const unsigned step
);

//! For convenience
unsigned ind_size(
    const unsigned beginning
,   const unsigned end
,   const unsigned step
);

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
);

#endif // UTILITIES_H_INCLUDED
