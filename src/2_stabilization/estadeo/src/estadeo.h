#ifndef ESTADEO_H
#define ESTADEO_H

#define LOCAL_MATRIX_BASED_SMOOTHING 0
#define LOCAL_LINEAR_MATRIX_BASED_SMOOTHING 1
#define LOCAL_LINEAR_POINT_BASED_SMOOTHING 2
#define N_SMOOTH_METHODS LOCAL_LINEAR_POINT_BASED_SMOOTHING


#include "utils.h"


/**
 *
 * Class for online video stabilization
 * Implements a circular array
 *
**/
class estadeo {

  public:
  
    estadeo(
      int   strat, //motion smoothing strategy
      int   np,    //number of parameters of the transformations
      float sigm,  //Gaussian standard deviation for smoothing
      int   verb   //switch on verbose mode
    );
    
    ~estadeo();
    
    void process_frame(
      float *I1,    //input previous image of video
      float *I2,    //input last image of video
      float *Ic,    //input last color image to warp
      Timer &timer, //keep runtimes
      int   nx,     //number of columns 
      int   ny,     //number of rows
      int   nz      //number of channels
    );
    
    float *get_H();

    float *get_smooth_H();
    
    int obtain_radius(){return (int)3*sigma;}

    
  private:
  
    void compute_motion(
      float *I1, //first image
      float *I2, //second image
      int   nx,  //number of columns
      int   ny   //number of rows
    );

    void motion_smoothing();
    
    void frame_warping(
      float *I, //frame to be warped
      int   nx, //number of columns   
      int   ny, //number of rows
      int   nz  //number of channels
    );

    
    //functions for motion smoothing
    void global_gaussian(int i);
    
    void local_matrix_based_smoothing();
    

  private:
  
    int   type;    //motion smoothing method
    int   Np;      //number of parameters in the transformation
    float sigma;   //Gaussian standard deviation 
    int   radius;  //radius of the Gaussian convolution
    float *Hs;     //last smoothing transform
    float *Hp;     //last stabilizing transform
    int   verbose; //swith on verbose mode
    
    //variables for the circular array
    int   N;    //circular array size
    int   Nf;   //number of frames
    int   fc;   //current frame position
    float *H;   //original matrix transformation
    float *Hc;  //composition of transformations
    float *H_1; //inverse transformations
};


#endif
