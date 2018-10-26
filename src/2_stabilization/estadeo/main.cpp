// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017-2018, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm> 

#include "estadeo.h"
#include "utils.h"
#include "transformation.h"

extern "C" {
#include "vpp.h"
}

#define PAR_DEFAULT_OUTVIDEO "output_video.raw"
#define PAR_DEFAULT_TRANSFORM SIMILARITY_TRANSFORM
#define PAR_DEFAULT_SMOOTHING LOCAL_MATRIX_BASED_SMOOTHING
#define PAR_DEFAULT_SIGMA_T 30.0
#define PAR_DEFAULT_OUTTRANSFORM "transform.mat"
#define PAR_DEFAULT_VERBOSE 0


/**
 *
 *  Print a help message 
 *
 */
void print_help(char *name)
{
  printf("\n  Usage: %s raw_input_video [OPTIONS] \n\n",
         name);
  printf("  Video stabilization:\n");
  printf("  'raw_input_video' is a video file in raw format (rgb24).\n");
  printf("  -----------------------------------------------\n");
  printf("  Converting to raw data:\n");
  printf("  'avconv -i video.mp4 -f rawvideo -pix_fmt rgb24 -y "
         "raw_video.raw'\n");
  printf("  to convert an mp4 video to raw format.\n");
  printf("  'avconv -f rawvideo -pix_fmt rgb24 -video_size 640x360 "
         "-framerate\n");
  printf("  30 -i output_video.raw -pix_fmt yuv420p output_video.mp4'\n");
  printf("  to convert a raw video to mp4 format.\n");
  printf("  -----------------------------------------------\n");
  printf("  More information in http://www.ipol.im \n\n");
  printf("  OPTIONS:\n"); 
  printf("  --------\n");
  printf("   -o name  output video name to write the computed raw video\n");
  printf("              default value '%s'\n", PAR_DEFAULT_OUTVIDEO);
  printf("   -t N     transformation type to be computed:\n");
  printf("              2.translation; 3.Euclidean transform;\n");
  printf("              4.similarity; 6.affinity; 8.homography\n"); 
  printf("              default value %d\n", PAR_DEFAULT_TRANSFORM);
  printf("   -m N     motion smoothing strategy:\n");
  printf("              0.local matrix-based smoothing;\n");
  printf("              1.local linear matrix-based smoothing\n");
  printf("              2.local linear point-based smoothing\n");
  printf("              default value %d\n", PAR_DEFAULT_SMOOTHING);
  printf("   -st N      Gaussian standard deviation for temporal dimension\n");
  printf("              default value %f\n", PAR_DEFAULT_SIGMA_T);
  printf("   -w name  write transformations to file\n");
  printf("   -f name  write stabilizing transformations to file\n");
  printf("   -v       switch on verbose mode \n\n\n");
}

/**
 *
 *  Read command line parameters 
 *
 */
int read_parameters(
                    int   argc, 
                    char  *argv[], 
                    char  **video_in,
                    char  *video_out,
                    char  **out_transform,
                    char  **out_smooth_transform,
                    int   &nparams,
                    int   &smooth_strategy,
                    float &sigma,
                    int   &verbose
                   )
{
  if (argc < 2){
    print_help(argv[0]); 
    return 0;
  }
  else{
    int i=1;
    *video_in=argv[i++];

    *out_transform=NULL;
    *out_smooth_transform=NULL;

    //assign default values to the parameters
    strcpy(video_out,PAR_DEFAULT_OUTVIDEO);
    nparams=PAR_DEFAULT_TRANSFORM;
    smooth_strategy=PAR_DEFAULT_SMOOTHING;
    sigma=PAR_DEFAULT_SIGMA_T;
    verbose=PAR_DEFAULT_VERBOSE;

    //read each parameter from the command line
    while(i<argc)
    {
      if(strcmp(argv[i],"-o")==0)
        if(i<argc-1)
          strcpy(video_out,argv[++i]);

      if(strcmp(argv[i],"-t")==0)
        if(i<argc-1)
          nparams=atof(argv[++i]);

      if(strcmp(argv[i],"-m")==0)
        if(i<argc-1)
          smooth_strategy=atoi(argv[++i]);

      if(strcmp(argv[i],"-st")==0)
        if(i<argc-1)
          sigma=atof(argv[++i]);

      if(strcmp(argv[i],"-w")==0)
        if(i<argc-1)
          *out_transform=argv[++i];

      if(strcmp(argv[i],"-f")==0)
        if(i<argc-1)
          *out_smooth_transform=argv[++i];

      if(strcmp(argv[i],"-v")==0)
        verbose=1;

      i++;
    }

    //check parameter values
    if(nparams!=2 && nparams!=3 && nparams!=4 && 
       nparams!=6 && nparams!=8) nparams=PAR_DEFAULT_TRANSFORM;
    if(smooth_strategy<0 || smooth_strategy>N_SMOOTH_METHODS)
      smooth_strategy=PAR_DEFAULT_SMOOTHING;
    if(sigma<0.01)
      sigma=0.01;
  }

  return 1;
}


/**
 *
 *  Function for converting an rgb image to grayscale levels
 * 
 **/
void rgb2gray(
              float *rgb,  //input color image
              float *gray, //output grayscale image
              int nx,      //number of pixels
              int ny, 
              int nz
             )
{
  int size=nx*ny;
  if(nz>=3)
#pragma omp parallel for
    for(int i=0;i<size;i++)
      gray[i]=(0.2989*rgb[i*nz]+0.5870*rgb[i*nz+1]+0.1140*rgb[i*nz+2]);
  else
#pragma omp parallel for
    for(int i=0;i<size;i++)
      gray[i]=rgb[i];
}


/**
 *
 *  Main program:
 *   This program reads the parameters from the console and
 *   then call the video stabilization method
 *
 */
int main (int argc, char *argv[])
{
  //parameters of the method
  char  *video_in, video_out[300];
  char  *out_transform, *out_stransform;
  int   width, height, nchannels=3;
  int   nparams, strategy, verbose;
  float sigma;

  //read the parameters from the console
  int result=read_parameters(
                             argc, argv, &video_in, video_out, &out_transform, &out_stransform,
                             nparams, strategy, sigma, verbose
                            );

  if(!result)
  {
    return EXIT_FAILURE;
  }

  if(verbose)
    printf(
           " Input video: '%s'\n Output video: '%s'\n"
           " Transformation: %d, Smoothing: %d\n",
           video_in, video_out, nparams, strategy
          );

  FILE* input = vpp_init_input(video_in, &width, &height, &nchannels);
  if (!input) {
    fprintf(stderr, "Error: Cannot initialize input '%s'.\n", video_in);
    return EXIT_FAILURE;
  }

  FILE* output = vpp_init_output(video_out, width, height, nchannels);
  if (!output) {
    fprintf(stderr, "Error: Cannot initialize output '%s'.\n", video_out);
    return EXIT_FAILURE;
  }

  int fsize=width*height;
  int csize=fsize*nchannels;

  //convert the input video to float and gray levels
  float *Ic=new float[csize];
  float *I1=new float[fsize];
  float *I2=new float[fsize];

  if(verbose) printf("\n Starting the stabilization\n");

  //read the first frame from input stream
  if (!vpp_read_frame(input, Ic, width, height, nchannels))
    return 1;
  if (!vpp_write_frame(output, Ic, width, height, nchannels))
    return 1;

  //convert it to grayscale
  rgb2gray(Ic, I1, width, height, nchannels);

  Timer timer;
  estadeo stab(strategy, nparams, sigma, verbose);

  int f = 0;
  while (true) {
    if (!vpp_read_frame(input, Ic, width, height, nchannels))
      break;

    //convert it to grayscale
    rgb2gray(Ic, I2, width, height, nchannels);

    //call the method for stabilizing current frame
    stab.process_frame(I1, I2, Ic, timer, width, height, nchannels);

    if(verbose) timer.print_time(f);

    //save the stabilized video to the output stream
    if (!vpp_write_frame(output, Ic, width, height, nchannels))
      break;

    std::swap(I1, I2);

    //save the motion transformations 
    // TODO: vpp
    if(out_transform!=NULL)
      save_transform(out_transform, stab.get_H(), nparams);

    //save the stabilizing transformation
    // TODO: vpp
    if(out_stransform!=NULL)
      save_transform(out_stransform, stab.get_smooth_H(), nparams);

    f++;
  }

  if(verbose) timer.print_avg_time(f);

  delete []Ic;
  delete []I1;
  delete []I2;

  return EXIT_SUCCESS;
}

