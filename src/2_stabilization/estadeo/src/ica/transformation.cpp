// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include "transformation.h"

#include <math.h>
#include <vector>


using namespace std;

/**
 *
 *  Function to compute the Jacobian matrix
 *  These parametrizations of the Jacobian are taken from the book of Zselinski
 *  (chapter 6 and 9)
 *
 */
void jacobian
(
  float *J,    //computed Jacobian
  vector<int> &p, //point coordinates
  int nparams, //number of parameters
  int nx       //number of columns
) 
{
  int N=p.size();
  switch(nparams) 
  {
    default: case TRANSLATION_TRANSFORM:  //p=(tx, ty) 
      for(int i=0; i<N; i++)
      {
        int c=2*i*nparams;
        J[c]  =1.0; J[c+1]=0.0;
        J[c+2]=0.0; J[c+3]=1.0;
      }
      break;
    case EUCLIDEAN_TRANSFORM:  //p=(tx, ty, tita)
      for(int i=0; i<N; i++)
      {
        int c=2*i*nparams;
        int x=p[i]%nx;
        int y=(int)(p[i]/nx);
        J[c]  =1.0; J[c+1]=0.0; J[c+2]=-y; 
        J[c+3]=0.0; J[c+4]=1.0; J[c+5]= x; 
      }
      break;
    case SIMILARITY_TRANSFORM: //p=(tx, ty, a, b)
      for(int i=0; i<N; i++)
      {
        int c=2*i*nparams;
        int x=p[i]%nx;
        int y=(int)(p[i]/nx); 
        J[c]  =1.0; J[c+1]=0.0; J[c+2]=x; J[c+3]=-y;
        J[c+4]=0.0; J[c+5]=1.0; J[c+6]=y; J[c+7]= x;
      }
      break;
    case AFFINITY_TRANSFORM:  //p=(tx, ty, a00, a01, a10, a11)
      for(int i=0; i<N; i++)
      {
        int c=2*i*nparams;
        int x=p[i]%nx;
        int y=(int)(p[i]/nx);
        J[c]  =1.0;J[c+1]=0.0;J[c+2]=  x;J[c+3]=  y;J[c+ 4]=0.0;J[c+ 5]=0.0;
        J[c+6]=0.0;J[c+7]=1.0;J[c+8]=0.0;J[c+9]=0.0;J[c+10]=  x;J[c+11]=  y;
      }
      break;    
    case HOMOGRAPHY_TRANSFORM: //p=(h00, h01,..., h21)
      for(int i=0; i<N; i++)
      {
        int c=2*i*nparams;
        int x=p[i]%nx;
        int y=(int)(p[i]/nx);
        J[c]  =  x; J[c+1]=  y; J[c+2]=1.0;  J[c+3]=0.0; 
        J[c+4]=0.0; J[c+5]=0.0; J[c+6]=-x*x; J[c+7]=-x*y;
        J[c+8]=0.0; J[c+9] =0.0; J[c+10]= 0.0; J[c+11]=   x; 
        J[c+12]= y; J[c+13]=1.0; J[c+14]=-x*y; J[c+15]=-y*y;
      }
      break;
  }
}

/**
 *
 *  Function to update the current transform with the computed increment
 *  x'(x;p) = x'(x;p) o x'(x;dp)^-1
 *
 */
void update_transform
(
  float *p,  //output accumulated transform
  float *dp, //computed increment
  int nparams //number of parameters
)
{
  switch(nparams) 
  {
    default: case TRANSLATION_TRANSFORM: //p=(tx, ty)
      for(int i = 0; i < nparams; i++) 
        p[i]-=dp[i];
      break;
    case EUCLIDEAN_TRANSFORM: //p=(tx, ty, tita)
    {
      float a=cos(dp[2]);
      float b=sin(dp[2]);
      float c=dp[0];
      float d=dp[1];
      float ap=cos(p[2]);
      float bp=sin(p[2]);
      float cp=p[0];
      float dp=p[1];
      float cost=a*ap+b*bp;
      float sint=a*bp-b*ap;
      p[0]=cp-bp*(b*c-a*d)-ap*(a*c+b*d);
      p[1]=dp-bp*(a*c+b*d)+ap*(b*c-a*d);
      p[2]=atan2(sint,cost);   
    }
    break;
    case SIMILARITY_TRANSFORM: //p=(tx, ty, a, b)
    {
      float a=dp[2];
      float b=dp[3];
      float c=dp[0];
      float d=dp[1];
      float det=(2*a+a*a+b*b+1);
      if(det*det>1E-10)
      {
        float ap=p[2];
        float bp=p[3];
        float cp=p[0];
        float dp=p[1];
        
        p[0]=cp-bp*(-d-a*d+b*c)/det+(ap+1)*(-c-a*c-b*d)/det;
        p[1]=dp+bp*(-c-a*c-b*d)/det+(ap+1)*(-d-a*d+b*c)/det;
        p[2]=b*bp/det+(a+1)*(ap+1)/det-1;
        p[3]=-b*(ap+1)/det+bp*(a+1)/det ;
      }
    }
    break;
    case AFFINITY_TRANSFORM: //p=(tx, ty, a00, a01, a10, a11)
    {
      float a=dp[2];
      float b=dp[3];
      float c=dp[0];
      float d=dp[4];
      float e=dp[5];
      float f=dp[1];
      float det=(a-b*d+e+a*e+1);
      if(det*det>1E-10)
      {
        float ap=p[2];
        float bp=p[3];
        float cp=p[0];
        float dp=p[4];
        float ep=p[5];
        float fp=p[1];
        
        p[0]=cp+(-f*bp-a*f*bp+c*d*bp)/det+(ap+1)*(-c+b*f-c*e)/det;
        p[1]=fp+dp*(-c+b*f-c*e)/det+(-f+c*d-a*f-f*ep-a*f*ep+d*d*ep)/det;
        p[2]=((1+ap)*(1+e)-d*bp)/det-1;
        p[3]=(bp+a*bp-b-b*ap)/det;
        p[4]=(dp*(1+e)-d-d*ep)/det;
        p[5]=(a+ep+a*ep+1-b*dp)/det-1; 
      }
    }
    break;
    case HOMOGRAPHY_TRANSFORM:   //p=(h00, h01,..., h21)
    {
      float a=dp[0];
      float b=dp[1];
      float c=dp[2];
      float d=dp[3];
      float e=dp[4];
      float f=dp[5];
      float g=dp[6];
      float h=dp[7];
      float ap=p[0];
      float bp=p[1];
      float cp=p[2];
      float dp=p[3];
      float ep=p[4];
      float fp=p[5];
      float gp=p[6];
      float hp=p[7];
        
      float det=f*hp+a*f*hp-c*d*hp+gp*(c-b*f+c*e)-a+b*d-e-a*e-1;
      if(det*det>1E-10)
      {
        p[0]=((d*bp-f*g*bp)+cp*(g-d*h+g*e)+(ap+1)*(f*h-e-1))/det-1;
        p[1]=(h*cp+a*h*cp-b*g*cp-bp-a*bp+c*g*bp+b-c*h+b*ap-c*h*ap)/det;
        p[2]=(f*bp+a*f*bp-c*d*bp+(ap+1)*(c-b*f+c*e)+cp*(-a+b*d-e-a*e-1))/det;
        p[3]=(fp*(g-d*h+g*e)+d-f*g+d*ep-f*g*ep+dp*(f*h-e-1))/det;
        p[4]=(b*dp-c*h*dp+h*fp+a*h*fp-b*g*fp-a+c*g-ep-a*ep+c*g*ep-1)/det-1;
        p[5]=(dp*(c-b*f+c*e)+f+a*f-c*d+f*ep+a*f*ep-c*d*ep+fp*(-a+b*d-e-a*e-1))/det;
        p[6]=(d*hp-f*g*hp+g-d*h+g*e+gp*(f*h-e-1))/det;
        p[7]=(h+a*h-b*g+b*gp-c*h*gp-hp-a*hp+c*g*hp)/det;          
      }
    }
    break;
  }
}


/**
 *
 *  Function to transform a 2D point (x,y) through a parametric model
 *
 */
void project
(
  int x,      //x component of the 2D point
  int y,      //y component of the 2D point
  float *p,  //parameters of the transformation
  float &xp, //x component of the transformed point
  float &yp, //y component of the transformed point
  int nparams //number of parameters
)
{
  switch(nparams) {
    default: case TRANSLATION_TRANSFORM: //p=(tx, ty) 
      xp=x+p[0];
      yp=y+p[1];
      break;
    case EUCLIDEAN_TRANSFORM:   //p=(tx, ty, tita)
      xp=cos(p[2])*x-sin(p[2])*y+p[0];
      yp=sin(p[2])*x+cos(p[2])*y+p[1];
      break;
    case SIMILARITY_TRANSFORM:  //p=(tx, ty, a, b)
      xp=(1+p[2])*x-p[3]*y+p[0];
      yp=p[3]*x+(1+p[2])*y+p[1];
      break;
    case AFFINITY_TRANSFORM:    //p=(tx, ty, a00, a01, a10, a11)
      xp=(1+p[2])*x+p[3]*y+p[0];
      yp=p[4]*x+(1+p[5])*y+p[1];
      break;
    case HOMOGRAPHY_TRANSFORM:  //p=(h00, h01,..., h21)
      float d=p[6]*x+p[7]*y+1;
      xp=((1+p[0])*x+p[1]*y+p[2])/d;
      yp=(p[3]*x+(1+p[4])*y+p[5])/d;
      break;
  }
}

/**
 *
 *  Function to convert a parametric model to its matrix representation
 *
 */
void params2matrix
(
  float *p,      //input parametric model
  float *matrix, //output matrix
  int   nparams  //number of parameters
)
{
  matrix[0]=matrix[4]=matrix[8]=1;
  matrix[1]=matrix[2]=matrix[3]=0;
  matrix[5]=matrix[6]=matrix[7]=0;
  switch(nparams) {
    default: case TRANSLATION_TRANSFORM: //p=(tx, ty) 
      matrix[2]=p[0];
      matrix[5]=p[1];
      break;
    case EUCLIDEAN_TRANSFORM:   //p=(tx, ty, tita)
      matrix[0]=cos(p[2]);
      matrix[1]=-sin(p[2]);
      matrix[2]=p[0];
      matrix[3]=sin(p[2]);
      matrix[4]=cos(p[2]);
      matrix[5]=p[1];
     break;
    case SIMILARITY_TRANSFORM:  //p=(tx, ty, a, b)
      matrix[0]=1+p[2];
      matrix[1]=-p[3];
      matrix[2]=p[0];
      matrix[3]=-matrix[1];
      matrix[4]=matrix[0];
      matrix[5]=p[1];
      break;
    case AFFINITY_TRANSFORM:    //p=(tx, ty, a00, a01, a10, a11)
      matrix[0]=1+p[2];
      matrix[1]=p[3];
      matrix[2]=p[0];
      matrix[3]=p[4];
      matrix[4]=1+p[5];
      matrix[5]=p[1];
      break;
    case HOMOGRAPHY_TRANSFORM:  //p=(h00, h01,..., h21)
      matrix[0]=1+p[0];
      matrix[1]=p[1];
      matrix[2]=p[2];
      matrix[3]=p[3];
      matrix[4]=1+p[4];
      matrix[5]=p[5];
      matrix[6]=p[6];
      matrix[7]=p[7];
      break;
    
  }
}


/**
 *
 *  Function to convert a matrix to its parametric model
 *
 */
void matrix2params
(
  float *matrix,
  float *p,
  int nparams
)
{
  float m[9],sint,cost;
  for(int i=0;i<9;i++) m[i]=matrix[i]/matrix[8];
  switch(nparams) {
    default: case TRANSLATION_TRANSFORM: //p=(tx, ty) 
      p[0]=m[2];
      p[1]=m[5];
      break;
    case EUCLIDEAN_TRANSFORM:   //p=(tx, ty, tita)
      p[0]=m[2];
      p[1]=m[5];
      cost=(m[0]+m[4])/2;   
      sint=(m[3]-m[1])/2;   
      p[2]=atan2(sint,cost);
      break;
    case SIMILARITY_TRANSFORM:  //p=(tx, ty, a, b)
      p[0]=m[2];
      p[1]=m[5];
      p[2]=(m[0]+m[4])/2-1;   
      p[3]=(m[3]-m[1])/2;   
      break;
    case AFFINITY_TRANSFORM:    //p=(tx, ty, a00, a01, a10, a11)
      p[0]=m[2];
      p[1]=m[5];
      p[2]=m[0]-1;
      p[3]=m[1];
      p[4]=m[3];
      p[5]=m[4]-1;
      break;
    case HOMOGRAPHY_TRANSFORM:  //p=(h00, h01,..., h21)
      p[0]=m[0]-1;
      p[1]=m[1];
      p[2]=m[2];
      p[3]=m[3];
      p[4]=m[4]-1;
      p[5]=m[5];
      p[6]=m[6];
      p[7]=m[7];
      break;
  }
}







/**
 *
 *  Function to compose two transforms 
 *  W(x;p) = W(x;p1) o W(x;p2)
 *
 */
void compose_transform
(
  float *p1,        //first input transform
  float *p2,        //second input transform
  float *p,         //output transform
  const int nparams  //number of parameters
)
{
  switch(nparams) 
  {
    default: case TRANSLATION_TRANSFORM: //p=(tx, ty)
      for(int i = 0; i < nparams; i++) 
	     p[i]=p1[i]+p2[i];
      break;
    case EUCLIDEAN_TRANSFORM:  //p=(tx, ty, tita)
    {
      const float a=cos(p1[2]);
      const float b=sin(p1[2]);
      const float c=p1[0];
      const float d=p1[1];
      const float ap=cos(p2[2]);
      const float bp=sin(p2[2]);
      const float cp=p2[0];
      const float dp=p2[1];
      const float cost=a*ap-b*bp;
      const float sint=a*bp+b*ap;
      p[0]=c+cp*a-dp*b;
      p[1]=d+dp*a+cp*b;
      p[2]=atan2(sint,cost);   
    }
    break;
    case SIMILARITY_TRANSFORM:  //p=(tx, ty, a, b)
    {
      const float a=p1[2];
      const float b=p1[3];
      const float c=p1[0];
      const float d=p1[1];
      const float ap=p2[2];
      const float bp=p2[3];
      const float cp=p2[0];
      const float dp=p2[1];
      p[0]=c-b*dp+cp*(a+1);
      p[1]=d+b*cp+dp*(a+1);
      p[2]=-b*bp+(a+1)*(ap+1)-1;
      p[3]=b*(ap+1)+bp*(a+1);
    }
    break;
    case AFFINITY_TRANSFORM:    //p=(tx, ty, a00, a01, a10, a11)
    {
      const float a=p1[2];
      const float b=p1[3];
      const float c=p1[0];
      const float d=p1[4];
      const float e=p1[5];
      const float f=p1[1];
      const float ap=p2[2];
      const float bp=p2[3];
      const float cp=p2[0];
      const float dp=p2[4];
      const float ep=p2[5];
      const float fp=p2[1];	
      p[0]=c+b*fp+cp*(a+1);
      p[1]=f+d*cp+fp*(e+1);
      p[2]=b*dp+(a+1)*(ap+1)-1;
      p[3]=b*(ep+1)+bp*(a+1);
      p[4]=d*(ap+1)+dp*(e+1);
      p[5]=d*bp+(ep+1)*(e+1)-1; 
    }
    break;
    case HOMOGRAPHY_TRANSFORM:   //p=(h00, h01,..., h21)
    {
      const float a=p1[0];
      const float b=p1[1];
      const float c=p1[2];
      const float d=p1[3];
      const float e=p1[4];
      const float f=p1[5];
      const float g=p1[6];
      const float h=p1[7];
      const float ap=p2[0];
      const float bp=p2[1];
      const float cp=p2[2];
      const float dp=p2[3];
      const float ep=p2[4];
      const float fp=p2[5];
      const float gp=p2[6];
      const float hp=p2[7];
      const float det=cp*g+fp*h+1; 
      
      if(det*det>1E-10)
      {
	     p[0]=(b*dp+c*gp+(a+1)*(ap+1))/det-1;
	     p[1]=(c*hp+b*(ep+1)+bp*(a+1))/det;
        p[2]=(c+a*cp+b*fp+cp)/det;
	     p[3]=(d*(ap+1)+f*gp+dp*(e+1))/det;
	     p[4]=(bp*d+f*hp+(ep+1)*(e+1))/det-1;
	     p[5]=(f+cp*d+fp*(e+1))/det;
	     p[6]=(gp+g*(ap+1)+dp*h)/det;
	     p[7]=(hp+h*(ep+1)+bp*g)/det;
      }
      else p[0]=p[1]=p[2]=p[3]=p[4]=p[5]=p[6]=p[7]=0;     
    }
    break;
   }
}



/**
 *
 *  Function to invert a transforms 
 *  W(x;p) = W(x;p1) o W(x;p2)
 *
 */
void inverse_transform
(
  float *p1,        //first input transform
  float *p,         //output transform
  const int nparams  //number of parameters
)
{
  switch(nparams) 
  {
    default: case TRANSLATION_TRANSFORM: //p=(tx, ty)
      p[0]=-p1[0];
      p[1]=-p1[1];
      break;
    case EUCLIDEAN_TRANSFORM:  //p=(tx, ty, tita)
    {
      const float a=p1[0];
      const float b=p1[1];
      const float c=p1[2];
      p[0]=-a*cos(c)-b*sin(c);
      p[1]= a*sin(c)-b*cos(c);
      p[2]=-c;
      break;
    }
    case SIMILARITY_TRANSFORM:  //p=(tx, ty, a, b)
    {
      const float a=p1[2];
      const float b=p1[3];
      const float c=p1[0];
      const float d=p1[1];
      const float det=2*a+a*a+b*b+1;
      if(det*det>1E-10)
      {
        p[0]=(-c-a*c-b*d)/det;
        p[1]=(-d-a*d+b*c)/det;
        p[2]=(a+1)/det-1;
        p[3]=-b/det;
      }
      else p[0]=p[1]=p[2]=p[3]=0;
    }
    break;
    case AFFINITY_TRANSFORM:    //p=(tx, ty, a00, a01, a10, a11)
    {
      const float a=p1[2];	
      const float b=p1[3];
      const float c=p1[0];
      const float d=p1[4];
      const float e=p1[5];
      const float f=p1[1];
      const float det=a-b*d+e+a*e+1;	
      if(det*det>1E-10)
      {
        p[0]=(-c+b*f-c*e)/det;
        p[1]=(-f-a*f+c*d)/det;
        p[2]=(e+1)/det-1;
        p[3]=-b/det;
        p[4]=-d/det;
        p[5]=(a+1)/det-1; 
      }
      else p[0]=p[1]=p[2]=p[3]=p[4]=p[5]=0;
    }
    break;
    case HOMOGRAPHY_TRANSFORM:   //p=(h00, h01,..., h21)
    {
      const float a=p1[0];
      const float b=p1[1];
      const float c=p1[2];
      const float d=p1[3];
      const float e=p1[4];
      const float f=p1[5];
      const float g=p1[6];
      const float h=p1[7];
      const float det=(-a+b*d-e-a*e-1);
      if(det*det>1E-10)
      {
        p[0]=(f*h-e-1)/det-1;
        p[1]=(b-c*h)/det;
        p[2]=(c-b*f+c*e)/det;
        p[3]=(d-f*g)/det;
        p[4]=(-a+c*g-1)/det-1;
        p[5]=(f+a*f-c*d)/det;
        p[6]=(g-d*h+g*e)/det;
        p[7]=(h+a*h-b*g)/det;
      }
      else p[0]=p[1]=p[2]=p[3]=p[4]=p[5]=p[6]=p[7]=0;
    }
    break;

  }
}





