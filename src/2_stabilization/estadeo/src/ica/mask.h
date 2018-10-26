// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López  <nmonzon@ctim.es>
// All rights reserved.

#ifndef MASK_H
#define MASK_H

#include <vector>

/**
 *
 * Compute the gradient with central differences
 *
 */
void gradient(
  float *input,  //input image
  float *dx,     //computed x derivative
  float *dy,     //computed y derivative
  int nx,        //image width
  int ny         //image height
);


/**
  *
  * Function to compute the gradient at given points
  * It does not treat border pixels
  *
**/
void gradient(
    float *input, //input image
    float *dx,    //computed x derivative
    float *dy,    //computed y derivative
    std::vector<int> &x, //positions to compute the gradient
    const int nx  //image width
);



/**
 *
 * Convolution with a Gaussian
 *
 */
void
gaussian (
  float *I,     //input image
  float *Is,    //output image
  int xdim,     //image width
  int ydim,     //image height
  float sigma,  //Gaussian sigma
  int precision = 4 //defines the size of the window
);

#endif
