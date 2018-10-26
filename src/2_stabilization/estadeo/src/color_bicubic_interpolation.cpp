// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015-2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// Copyright (C) 2014, Nelson Monzón López <nmonzon@ctim.es>
// All rights reserved.


#include "color_bicubic_interpolation.h"
#include "bicubic_interpolation.h"
#include "transformation.h"


/**
  *
  * Neumann boundary condition test
  *
**/
int neumann_bc (int x, int nx)
{
  if (x<0) x=0;
  else if (x >= nx)
    x = nx - 1;
  return x;
}


/**
  *
  * Compute the bicubic interpolation of a point in an image. 
  * It detects if the point goes beyond the image domain
  *
**/
float bicubic_interpolation(
  float *input, //image to be interpolated
  float uu,     //x component of the vector field
  float vv,     //y component of the vector field
  int nx,       //width of the image
  int ny,       //height of the image
  int nz,       //number of channels of the image
  int k         //actual channel
)
{
  if(uu>nx || uu<-1 || vv>ny || vv<-1) 
    return 0;
  else
  {
     int sx = (uu < 0) ? -1 : 1;
     int sy = (vv < 0) ? -1 : 1;

     int x, y, mx, my, dx, dy, ddx, ddy;
  
     x = neumann_bc ((int) uu, nx);
     y = neumann_bc ((int) vv, ny);
     mx = neumann_bc ((int) uu - sx, nx);
     my = neumann_bc ((int) vv - sy, ny);
     dx = neumann_bc ((int) uu + sx, nx);
     dy = neumann_bc ((int) vv + sy, ny);
     ddx = neumann_bc ((int) uu + 2 * sx, nx);
     ddy = neumann_bc ((int) vv + 2 * sy, ny);

     //obtain the interpolation points of the image
     float p11 = input[(mx  + nx * my) * nz + k];
     float p12 = input[(x   + nx * my) * nz + k];
     float p13 = input[(dx  + nx * my) * nz + k];
     float p14 = input[(ddx + nx * my) * nz + k];

     float p21 = input[(mx  + nx * y) * nz + k];
     float p22 = input[(x   + nx * y) * nz + k];
     float p23 = input[(dx  + nx * y) * nz + k];
     float p24 = input[(ddx + nx * y) * nz + k];
     
     float p31 = input[(mx  + nx * dy) * nz + k];
     float p32 = input[(x   + nx * dy) * nz + k];
     float p33 = input[(dx  + nx * dy) * nz + k];
     float p34 = input[(ddx + nx * dy) * nz + k];

     float p41 = input[(mx  + nx * ddy) * nz + k];
     float p42 = input[(x   + nx * ddy) * nz + k];
     float p43 = input[(dx  + nx * ddy) * nz + k];
     float p44 = input[(ddx + nx * ddy) * nz + k];
     
     //create array
     float pol[4][4] = { 
       {p11, p21, p31, p41}, {p12, p22, p32, p42},
       {p13, p23, p33, p43}, {p14, p24, p34, p44}
     };

     //return interpolation
     return bicubic_interpolation(pol, (float) uu-x, (float) vv-y);
  }
}



/**
  *
  * Compute the bicubic interpolation of an image from a parametric trasform
  *
**/
void bicubic_interpolation(
  float *input,   //image to be warped
  float *output,  //warped output image with bicubic interpolation
  float *params,  //x component of the vector field
  int nparams,    //number of parameters of the transform
  int nx,         //width of the image
  int ny,         //height of the image
  int nz          //number of channels of the image       
)
{
  #pragma omp parallel for
  for (int i=0; i<ny; i++)
    for (int j=0; j<nx; j++)
    {
      int p=i*nx+j;
      float x, y;

      //transform coordinates using the parametric model
      project(j, i, params, x, y, nparams);
      
      //obtain the bicubic interpolation at position (uu, vv)
      for(int k=0; k<nz; k++)
        output[p*nz+k]=bicubic_interpolation(
          input, x, y, nx, ny, nz, k
        );
    }
}


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
)
{
  #pragma omp parallel for
  for(int i = 0; i < ny; i++)
    for(int j = 0; j < nx; j++)
    {
       float uu, vv;

       //transform coordinates using the parametric model
       project(j, i, params, uu, vv, nparams);
 
       if(uu<1 || uu>nx-2 || vv<1 || vv>ny-2)
         for(int k=0; k<nz; k++)
           output[(j+nx*i)*nz+k]=0;
       else {
         int sx=(uu<0)? -1: 1;
         int sy=(vv<0)? -1: 1;
         int x, y, dx, dy;

         x =(int) uu;
         y =(int) vv;
         dx=(int) uu+sx;
         dy=(int) vv+sy;
         
         for(int k=0; k<nz; k++){
           float p1=input[(x +nx*y)*nz+k];
           float p2=input[(dx+nx*y)*nz+k];
           float p3=input[(x +nx*dy)*nz+k];
           float p4=input[(dx+nx*dy)*nz+k];

           float e1=((float) sx*(uu-x));
           float E1=((float) 1.0-e1);
           float e2=((float) sy*(vv-y));
           float E2=((float) 1.0-e2);

           float w1=E1*p1+e1*p2;
           float w2=E1*p3+e1*p4;

           output[(j+nx*i)*nz+k]=E2*w1+e2*w2;
         }
       }
    }
}


