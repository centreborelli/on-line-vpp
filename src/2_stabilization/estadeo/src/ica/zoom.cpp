// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015-2018, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.


#include <math.h>
#include <stdio.h>

#include "zoom.h"
#include "mask.h"
#include "transformation.h"

#define ZOOM_SIGMA_ZERO 0.7


/**
  *
  * Compute the size of a zoomed image from the zoom factor
  *
**/
void zoom_size 
(
  int nx,   //width of the orignal image
  int ny,   //height of the orignal image          
  int &nxx, //width of the zoomed image
  int &nyy  //height of the zoomed image
)
{
  nxx = (int) ((float) nx/2.+0.5);
  nyy = (int) ((float) ny/2.+0.5);
}


/**
  *
  * Function to downsample the image
  *
**/
void zoom_out
(
  float *I,    //input image
  float *Iout, //output image
  int   nx,    //image width
  int   ny     //image height          
)
{
  int nxx, nyy; 
  float *Is=new float[nx*ny];

  //calculate the size of the zoomed image
  zoom_size(nx, ny, nxx, nyy);
  
  //compute the Gaussian sigma for smoothing
  //sigma=sqrt(1.4^2-0.7^2) //1.21
  float sigma=ZOOM_SIGMA_ZERO*sqrt(3); 

  //pre-smooth the image
  gaussian(I, Is, nx, ny, sigma, 4);
  
  //re-sample image
  #pragma omp parallel for
  for (int i1=0; i1<nyy; i1++)
    for (int j1=0; j1<nxx; j1++)
    {
      int i2=2*i1;
      int j2=2*j1;
      Iout[i1*nxx+j1]=Is[i2*nx+j2];
    }   
  delete []Is;
}


/**
  *
  * Function to upsample the parameters of the transformation
  *
**/
void zoom_in_parameters 
(
  float *p,    //input image
  float *pout, //output image   
  int nparams, //number of parameters
  int   nx,    //width of the original image
  int   ny,    //height of the original image
  int   nxx,   //width of the zoomed image
  int   nyy    //height of the zoomed image
)
{
  //compute the zoom factor
  float factorx=((float)nxx/nx);
  float factory=((float)nyy/ny);
  float nu=(factorx>factory)?factorx:factory;

  switch(nparams) {
    default: case TRANSLATION_TRANSFORM: //p=(tx, ty) 
      pout[0]=p[0]*nu;
      pout[1]=p[1]*nu;
      break;
    case EUCLIDEAN_TRANSFORM: //p=(tx, ty, tita)
      pout[0]=p[0]*nu;
      pout[1]=p[1]*nu;
      pout[2]=p[2];
      break;
    case SIMILARITY_TRANSFORM: //p=(tx, ty, a, b)
      pout[0]=p[0]*nu;
      pout[1]=p[1]*nu;
      pout[2]=p[2];
      pout[3]=p[3];
      break;
    case AFFINITY_TRANSFORM: //p=(tx, ty, a00, a01, a10, a11)
      pout[0]=p[0]*nu;
      pout[1]=p[1]*nu;
      pout[2]=p[2];
      pout[3]=p[3];
      pout[4]=p[4];
      pout[5]=p[5];
      break;
    case HOMOGRAPHY_TRANSFORM: //p=(h00, h01,..., h21)
      pout[0]=p[0];
      pout[1]=p[1];
      pout[2]=p[2]*nu;
      pout[3]=p[3];
      pout[4]=p[4];
      pout[5]=p[5]*nu;
      pout[6]=p[6]/nu;
      pout[7]=p[7]/nu;
      break;
  }
}
