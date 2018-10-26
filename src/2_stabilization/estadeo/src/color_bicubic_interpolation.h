// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.


#ifndef COLOR_BICUBIC_INTERPOLATION_H
#define COLOR_BICUBIC_INTERPOLATION_H


/**
  *
  * Compute the bicubic interpolation of a point in an image. 
  * Detects if the point goes outside the image domain
  *
**/
float
bicubic_interpolation(
  float *input,//image to be interpolated
  float uu,    //x component of the vector field
  float vv,    //y component of the vector field
  int nx,      //width of the image
  int ny,      //height of the image
  int nz,      //number of channels of the image
  int k        //actual channel
);


/**
  *
  * Compute the bicubic interpolation of an image from a parametric trasform
  *
**/
void bicubic_interpolation(
  float *input,        //image to be warped
  float *output,       //warped output image with bicubic interpolation
  float *params,       //x component of the vector field
  int nparams,         //number of parameters of the transform
  int nx,              //width of the image
  int ny,              //height of the image
  int nz               //number of channels of the image       
);


/**
  *
  * Function to warp the image using bilinear interpolation
  *
**/
void bilinear_interpolation(
  float *input,   //image to be warped
  float *output,  //warped output image with bicubic interpolation
  float *params,  //x component of the vector field
  int nparams,    //number of parameters of the transform
  int nx,         //width of the image
  int ny,         //height of the image
  int nz          //number of channels of the image       
);


#endif
