// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017-2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#include "estadeo.h"
#include "color_bicubic_interpolation.h"
#include "inverse_compositional_algorithm.h"
//#include "gaussian_conv_dct.h"
#include "transformation.h"
#include "matrix.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


estadeo::estadeo(int strategy, int np, float sigm, int verb): 
         type(strategy), Np(np), sigma(sigm), verbose(verb)
{
  radius=obtain_radius();
  N=(2*radius+1);
  Nf=1;
  fc=0;
 
  //allocate motion transformations for the circular array
  H  =new float[N*Np]; //motion transformations
  Hc =new float[N*Np]; //composition of transformations
  H_1=new float[N*Np]; //inverse transformations

  //introduce identity matrix for the first transform and its inverse
  for(int i=0; i<Np; i++) H[i]=H_1[i]=0;
  
  //allocate last smooth transform
  Hs=new float[Np];
  
  //allocate last stabilizing transform
  Hp=new float[Np];
}

estadeo::~estadeo()
{
  delete []H;
  delete []Hc;
  delete []H_1;
  delete []Hs;
  delete []Hp;
}

/**
  *
  * Main function for Online Video Estabilization
  * Process a frame each time
  * 
**/
void estadeo::process_frame(
  float *I1,    //input previous image of video
  float *I2,    //input last image of video
  float *Ic,    //input last color image to warp
  Timer &timer, //keep runtimes
  int   nx,     //number of columns 
  int   ny,     //number of rows
  int   nz      //number of channels
)
{ 
  //increase the number of frames
  Nf++;
  
  //increase the position of the current frame in the circular array
  fc++;
  if(fc>N-1) fc=0;
  
  //step 1. Compute motion between the last two frames
  if(verbose) timer.set_t1();
  compute_motion(I1, I2, nx, ny);
    
  //step 2. Smooth until the last transformation 
  if(verbose) timer.set_t2();
  motion_smoothing();
    
  //step 3. Warp the image
  if(verbose) timer.set_t3();
  frame_warping(Ic, nx, ny, nz);
  if(verbose) timer.set_t4();
}


/**
  *
  * Function for estimating the transformation between two frames
  *
**/
void estadeo::compute_motion
(
  float *I1, //first image
  float *I2, //second image
  int   nx,  //number of columns
  int   ny   //number of rows
)
{
  //parameters for the direct method
  int   nscales=100;
  float TOL=1E-3;
  float lambda=0;
  float N;
  int   robust=LORENTZIAN; //QUADRATIC; //
  int   max_d=200;
  int   min_d=50;

  N=1+log(((nx<ny)?nx:ny)/min_d)/log(2.);
  if ((int) N<nscales) nscales=(int) N;

  //motion estimation through direct methods
  pyramidal_inverse_compositional_algorithm(
    I1, I2, get_H(), Np, nx, ny, nscales, TOL, robust, lambda, max_d
  );
}


/**
  *
  * Function for online motion_smoothing
  * 
  *
**/
void estadeo::motion_smoothing()
{
  switch(type) 
  {
    default: case LOCAL_MATRIX_BASED_SMOOTHING:
      local_matrix_based_smoothing();
      break;
    /*case LOCAL_LINEAR_MATRIX_BASED_SMOOTHING:
      local_linear_matrix_based_smoothing(H, Hp, nparams, ntransforms, sigma);
      break;
    case LOCAL_LINEAR_POINT_BASED_SMOOTHING:
      local_linear_point_based_smoothing(H, Hp, nparams, ntransforms, sigma);
      break;*/
  }
}


/**
  *
  * Function for online warping the last frame of the video
  *
**/
void estadeo::frame_warping
(
  float *I, //frame to be warped
  int   nx, //number of columns   
  int   ny, //number of rows
  int   nz  //number of channels
)
{
  int size=nx*ny*nz;

  float *I2=new float[nx*ny*nz];

  //warp the image
  bicubic_interpolation(I, I2, Hp, Np, nx, ny, nz);
  //bilinear_interpolation(I, I2, Hp, Np, nx, ny, nz);

  //copy warped image
  for(int j=0; j<size; j++)
    I[j]=I2[j];

  delete []I2;
}



/**
  *
  * Function to return the motion of the last transformation
  *
**/
float *estadeo::get_H()
{
  return &H[fc*Np];
}


/**
  *
  * Function to compute the motion of the last smooth transformation
  *
**/
float *estadeo::get_smooth_H()
{
  float *H_1=new float[Np];
  float *Htmp=new float[Np];

  inverse_transform(Hp, H_1, Np);
  compose_transform(get_H(), Hp, Htmp, Np);
  compose_transform(H_1, Htmp, Hs, Np);
  
  delete []H_1;
  delete []Htmp;
  
  return Hs;
}


//Gaussian convolution
void estadeo::global_gaussian(int i)
{
  //obtain current radius
  int rad=radius; 
  if(rad>=Nf) rad=Nf-1;

  //obtain current size of circular array  
  int n=N;   
  if(n>Nf) n=Nf;

  //Gaussian convolution in  each parameter separately
  for(int p=0; p<Np; p++)
  {
    double average=0.0;
    double sum=0.0;
    
    for(int j=i-rad;j<=i+rad;j++)
    {
      double value=0;
      
      //Neumann boundary conditions
      if(j<0)
        value=Hc[-j*Np+p];
      else if(j>=Nf){
        int l=(2*Nf-1-j)%n;
        value=Hc[l*Np+p];
      }
      else
        value=Hc[(j%n)*Np+p];
      
      //increase accumulator
      double norm=0.5*(j-i)*(j-i)/(sigma*sigma);
      double gauss=exp(-norm);
      average+=gauss*value;
      sum+=gauss;
    }
    Hs[p]=(float) (average/sum);
  }
}


//local matrix based smoothing approach
//SE PUEDEN CALCULAR DE FORMA INCREMENTAL Hc (con una sola H_1)
void estadeo::local_matrix_based_smoothing()
{
  //obtain current radius
  int rad=radius; 
  if(rad>=Nf) rad=Nf-1;

  //obtain current size of circular array  
  int n=N;   
  if(n>Nf) n=Nf;

  //compute inverse transform
  inverse_transform(&(H[fc*Np]), &(H_1[fc*Np]), Np);

  //recompute the stabilization for past frames using the circular array
  for(int i=Nf-rad; i<Nf; i++)
  {
    int t1=(i-rad>0)?i-rad: 0;

    //compute backward transformations
    int f=i%n;     //current frame in circular array
    int l=(i-1)%n; //previous frame in circular array
    
    for(int j=0;j<Np;j++) 
      Hc[l*Np+j]=H_1[f*Np+j];
    for(int j=i-2;j>=t1;j--)
    {
      int l1=(j+1)%n;
      int l2=j%n;
      compose_transform(&(H_1[l1*Np]), &(Hc[l1*Np]), &(Hc[l2*Np]), Np);
    }

    //introduce the identity matrix in the current frame
    for(int j=0;j<Np;j++) Hc[f*Np+j]=0;

    //compute forward transformations
    if(i<Nf-1)
    { 
      int r=(i+1)%n;
      for(int j=0;j<Np;j++) 
        Hc[r*Np+j]=H[r*Np+j];
      for(int j=i+2;j<=Nf;j++)
      {
        int r1=j%n;
        int r2=(j-1)%n;
        compose_transform(&(H[r1*Np]), &(Hc[r2*Np]), &(Hc[r1*Np]), Np);  
      }
    }

    //smooth transforms with a discrete Gaussian kernel
    global_gaussian(i);

    //compute inverse transformations 
    inverse_transform(Hs, Hp, Np);
  }
}



/*//Gaussian convolution for local methods
void estadeo::local_gaussian
(
  float *H,          //original matrix transformations
  float *Hs,         //smooth output matrix transformations
  int   i,           //frame number
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video  
  float sigma,       //Gaussian standard deviation
  float *H_1         //inverse transforms
)
{
  int radius=obtain_radius(sigma);

  if(radius>=ntransforms) 
    radius=ntransforms-1;
  
  //Gaussian convolution in each parameter separately
  for(int p=0;p<nparams;p++)
  { 
    float average=0.0;
    float sum=0.0;
    
    for(int j=i-radius;j<=i+radius;j++)
    {
      float value=0;
      
      //test boundary conditions
      if(j<0)
        value=H_1[(-j)*nparams+p];
      else if(j>=ntransforms) 
        value=H_1[(2*ntransforms-1-j)*nparams+p];
      else 
        value=H[j*nparams+p];
      
      //increase accumulator
      float norm=0.5*(j-i)*(j-i)/(sigma*sigma);
      float gauss=exp(-norm);
      average+=gauss*value;
      sum+=gauss;
    }
    Hs[p]=average/sum;
  }
}

//Gaussian convolution with a set of points
float estadeo::point_gaussian(
  float *x,          //set of points
  int   i,           //frame number
  int   ntransforms, //number of transforms
  float sigma        //Gaussian standard deviation
)
{
  float average=0.0;
  float sum=0.0;
  
  int radius=obtain_radius(sigma);

  if(radius>=ntransforms) 
    radius=ntransforms-1;
  
  //Gaussian convolution
  for(int j=i-radius;j<=i+radius;j++)
  {
    float value=0;
    
    //test boundary conditions
    if(j<0)
      value=x[-j];
    else if(j>=ntransforms) 
      value=x[2*ntransforms-1-j];
    else 
      value=x[j];

    float dx=j-i;
    float norm=0.5*dx*dx/(sigma*sigma);
    float gauss=exp(-norm);
    average+=gauss*value;
    sum+=gauss;
  }
  return average/sum;
}


//Matrix DCT Gaussian convolution
void estadeo::matrix_gaussian_dct
(
  float *H,      //original matrix transformations
  float *Hs,     //smooth output matrix transformations
  int   nparams, //type of matrix transformation
  int   N,       //number of frames of the video  
  float sigma    //Gaussian standard deviation
)
{  
  num *dest=new num[3*N-2];
  num *src=new num[3*N-2];
  dct_coeffs c;
  
  //convolution in each matrix position
  for(int p=0; p<nparams; p++)
  {
    //copy the original image
    for(int i=N-1; i<2*N-1; i++)
      src[i]=H[(i-N+1)*nparams+p];
    
    //Neumann boundary conditions
    for(int i=0; i<N-1; i++)
      src[i]=H[(N-i-2)*nparams+p];
    for(int i=2*N-1; i<3*N-2; i++)
      src[i]=H[(3*N-i-2)*nparams+p];

    //apply DCT Gaussian convolution
    if (!(dct_precomp(&c, dest, src, 3*N-2, 1, sigma)))
      printf("Error in Gaussian convolution with DCT.");
    else {
      dct_gaussian_conv(c);
      dct_free(&c);
    }

    //copy the signal in the domain
    for(int i=0; i<N; i++)
      Hs[i*nparams+p]=dest[N-1+i];
  }

  delete []src;
  delete []dest;
}


//local linear matrix-based smoothing
void estadeo::local_linear_matrix_based_smoothing
(
  float *H,          //original matrix transformations
  float *Hp,         //smooth output matrix transformations
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video  
  float sigma        //Gaussian standard deviations
)
{
  float *HH=new float[ntransforms*9];
  float *Hi=new float[ntransforms*9];
  float *Hs=new float[ntransforms*9];

  //convert from params to matrices
  for(int i=1;i<ntransforms;i++) 
    params2matrix(&(H[i*nparams]), &(HH[i*9]), nparams);

  //identity matrix in the first position
  Hi[1]=Hi[2]=Hi[3]=0;
  Hi[5]=Hi[6]=Hi[7]=0;
  Hi[0]=Hi[4]=Hi[8]=1;

  //compute the virtual matrix trajectories
  for(int i=1;i<ntransforms;i++) 
  {
    for(int j=0;j<9;j++)
      //accumulate the homographies
      Hi[i*9+j]=Hi[(i-1)*9+j]+HH[i*9+j];
    
    //subtract the identity matrix
    Hi[i*9]-=1;Hi[i*9+4]-=1;Hi[i*9+8]-=1;
  }

  //convolve the virtual trajectories with a Gaussian kernel
  matrix_gaussian_dct(Hi, Hs, 9, ntransforms, sigma);

  for(int i=0;i<ntransforms;i++) 
  {
    //compute the correction matrix
    for(int j=0;j<9;j++)
      Hs[i*9+j]-=Hi[i*9+j];

    //add the identity matrix
    Hs[i*9]+=1;Hs[i*9+4]+=1;Hs[i*9+8]+=1;

    //convert homographies to params
    matrix2params(&(Hs[i*9]),&(Hp[i*nparams]),nparams);

    //compute its inverse 
    inverse_transform(&(Hp[i*nparams]), &(Hp[i*nparams]), nparams);        
  }
  
  delete []HH;
  delete []Hi;
  delete []Hs;
}


//Point DCT Gaussian convolution
void estadeo::point_gaussian_dct
(
  float *x,    //input set of points
  float *xs,   //output set of smoothed points 
  int   N,     //number of points
  float sigma  //Gaussian standard deviation
)
{  
  num *dest=new num[3*N-2];
  num *src=new num[3*N-2];
  dct_coeffs c;

  //copy the original signal
  for(int i=N-1; i<2*N-1; i++)
    src[i]=x[i-N+1];
  
  //Neumann boundary conditions
  for(int i=0; i<N-1; i++)
    src[i]=x[N-i-2];
  for(int i=2*N-1; i<3*N-2; i++)
    src[i]=x[3*N-i-2];

  //apply DCT Gaussian convolution
  if (!(dct_precomp(&c, dest, src, 3*N-2, 1, sigma)))
    printf("Error in Gaussian convolution with DCT.");
  else {
    dct_gaussian_conv(c);
    dct_free(&c);
  }

  //copy the signal in the domain
  for(int i=0; i<N; i++)
    xs[i]=dest[N-1+i];
  
  delete []src;
  delete []dest;
}


//local linear point based smoothing approach
void estadeo::local_linear_point_based_smoothing
(
  float *H,          //original matrix transformations
  float *Hp,         //smooth output matrix transformations
  int   nparams,     //type of matrix transformation
  int   ntransforms, //number of frames of the video  
  float sigma        //Gaussian standard deviations
)
{ 
  float *x0=new  float[ntransforms];
  float *x1=new  float[ntransforms];
  float *x2=new  float[ntransforms];
  float *x3=new  float[ntransforms];
  float *y0=new  float[ntransforms];
  float *y1=new  float[ntransforms];
  float *y2=new  float[ntransforms];
  float *y3=new  float[ntransforms];
  float *x0s=new float[ntransforms];
  float *x1s=new float[ntransforms];
  float *x2s=new float[ntransforms];
  float *x3s=new float[ntransforms];
  float *y0s=new float[ntransforms];
  float *y1s=new float[ntransforms];
  float *y2s=new float[ntransforms];
  float *y3s=new float[ntransforms];
  float *HH=new  float[ntransforms*9];

  //convert from params to matrices
  for(int i=0;i<ntransforms;i++) 
    params2matrix(&(H[i*nparams]), &(HH[i*9]), nparams);

  //choose four fixed points for all frames
  float xp[4]={0, 0, 500, 500};
  float yp[4]={0, 500, 0, 500};

  //tracking a set of points to be smoothed
  x0[0]=xp[0]; y0[0]=yp[0];
  x1[0]=xp[1]; y1[0]=yp[1];
  x2[0]=xp[2]; y2[0]=yp[2];
  x3[0]=xp[3]; y3[0]=yp[3];
  
  //compute the virtual trajectories
  for(int i=1;i<ntransforms;i++) 
  { 
    float dx, dy;
    Hx(&(HH[i*9]),xp[0],yp[0],dx,dy);
    x0[i]=x0[i-1]+(dx-xp[0]);
    y0[i]=y0[i-1]+(dy-yp[0]);

    Hx(&(HH[i*9]),xp[1],yp[1],dx,dy);
    x1[i]=x1[i-1]+(dx-xp[1]);
    y1[i]=y1[i-1]+(dy-yp[1]);

    Hx(&(HH[i*9]),xp[2],yp[2],dx,dy);
    x2[i]=x2[i-1]+(dx-xp[2]);
    y2[i]=y2[i-1]+(dy-yp[2]);

    Hx(&(HH[i*9]),xp[3],yp[3],dx,dy);
    x3[i]=x3[i-1]+(dx-xp[3]);
    y3[i]=y3[i-1]+(dy-yp[3]);
  }
  
  //DCT Gaussian convolution of each virtual trajectory
  point_gaussian_dct(x0, x0s, ntransforms, sigma);
  point_gaussian_dct(x1, x1s, ntransforms, sigma);
  point_gaussian_dct(x2, x2s, ntransforms, sigma);
  point_gaussian_dct(x3, x3s, ntransforms, sigma);
  point_gaussian_dct(y0, y0s, ntransforms, sigma);
  point_gaussian_dct(y1, y1s, ntransforms, sigma);
  point_gaussian_dct(y2, y2s, ntransforms, sigma);
  point_gaussian_dct(y3, y3s, ntransforms, sigma);
  
  for(int i=0;i<ntransforms;i++) 
  {
    x0s[i]+=xp[0]-x0[i];
    x1s[i]+=xp[1]-x1[i];
    x2s[i]+=xp[2]-x2[i];
    x3s[i]+=xp[3]-x3[i];
    y0s[i]+=yp[0]-y0[i];
    y1s[i]+=yp[1]-y1[i];
    y2s[i]+=yp[2]-y2[i];
    y3s[i]+=yp[3]-y3[i];
    
    //calculate the smoothed homography    
    float tmp[9];
    compute_H(
      x0s[i],x1s[i],x2s[i],x3s[i],
      y0s[i],y1s[i],y2s[i],y3s[i],      
      xp[0],xp[1],xp[2],xp[3],      
      yp[0],yp[1],yp[2],yp[3],     
      tmp
    );

    float Hout[9];
    for(int j=0; j<9; j++) Hout[j]=tmp[j];
    
    //convert homographies to params
    matrix2params(Hout,&(Hp[i*nparams]),nparams);
  }

  delete []x0s;
  delete []x1s;
  delete []x2s;
  delete []x3s;
  delete []y0s;
  delete []y1s;
  delete []y2s;
  delete []y3s;
  delete []x0;
  delete []x1;
  delete []x2;
  delete []x3;
  delete []y0;
  delete []y1;
  delete []y2;
  delete []y3;
  delete []HH;
}
*/


