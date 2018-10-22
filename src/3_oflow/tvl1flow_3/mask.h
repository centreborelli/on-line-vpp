// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.


#ifndef MASK_H
#define MASK_H

/**
 *
 * Details on how to compute the divergence and the grad(u) can be found in:
 * [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
 * Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
 *
 **/


/**
 *
 * Function to compute the divergence with backward differences
 * (see [2] for details)
 *
 **/
void divergence(
		const float *v1, // x component of the vector field
		const float *v2, // y component of the vector field
		float *div,      // output divergence
		const int nx,    // image width
		const int ny     // image height
	       );


/**
 *
 * Function to compute the gradient with forward differences
 * (see [2] for details)
 *
 **/
void forward_gradient(
		const float *f, //input image
		float *fx,      //computed x derivative
		float *fy,      //computed y derivative
		const int nx,   //image width
		const int ny    //image height
		);


/**
 *
 * Function to compute the gradient with centered differences
 *
 **/
void centered_gradient(
		const float *input,  //input image
		float *dx,           //computed x derivative
		float *dy,           //computed y derivative
		const int nx,        //image width
		const int ny         //image height
		);


/**
 *
 * In-place Gaussian smoothing of an image
 *
 */
void gaussian(
	float *I,             // input/output image
	const int xdim,       // image width
	const int ydim,       // image height
	const double sigma    // Gaussian sigma
);


#endif//MASK_H
