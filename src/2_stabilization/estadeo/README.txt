Estadeo Video Stabilization
***************************

*******
SUMMARY
*******

This program implements a video stabilization method using parametric models.
Motion smoothing is one of the key steps for video stabilization. Its objective
is to remove the undesired motion of a video from the transformations computed
in a previous motion estimation process. In this work, we study several 
smoothing strategies. We assume that the motion of the scene is captured
through a set of parametric models, which are the most commonly used 
techniques. We propose a classification of strategies based on convolutions with
a Gaussian kernel. These strategies are classified into global and local 
methods, depending on how the reference frame is defined. Our study follows a 
unified formal approach and explains some details for the implementation of the 
techniques, such as the definition of correct boundary conditions in each case. 
The experiments show that local methods usually outperform global ones but, 
under certain conditions, global methods have a similar performance.

This program is part of an IPOL publication (http://www.ipol.im/).

Reference articles:

[1] Javier Sanchez and Jean-Michel Morel, "Comparison of Motion Smoothing 
Strategies for Video Stabilization using Parametric Models", Image Processing
Online, 2017.

[2] Javier Sanchez and Jean-Michel Morel, "MOTION SMOOTHING STRATEGIES FOR 2D
VIDEO STABILIZATION", SIAM J. Imaging Sciences, 2017.


******
AUTHOR
******

Javier Sánchez Pérez <jsanchez@ulpgc.es> 
Centro de Tecnologías de la Imagen (CTIM) 
Universidad de Las Palmas de Gran Canaria


*******
VERSION
*******

Version 1, released on November 2, 2017


*******
LICENSE
*******

This program is free software: you can use, modify and/or redistribute it
under the terms of the simplified BSD License. You should have received a
copy of this license along this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.

Copyright (C) 2017, Javier Sánchez Pérez <jsanchez@ulpgc.es>
All rights reserved.


***********
COMPILATION
***********

Required environment: Any unix-like system with a standard compilation
environment (make and C and C++ compilers)

Compilation instructions: run "make" to produce two executables:
 - "estadeo" the main algorithm
 - "generate_output" auxiliary program for the online demo

 
*****
USAGE
*****

The program reads an input video in raw format and produces an output raw
video with the video stabilization. The input videos need to be converted
to raw format. This can be done using the 'avconv' or 'ffmpeg' programs as:

  1. Converting an mpeg 4 video file to raw data:

    'avconv -i video.mp4 -f rawvideo -pix_fmt rgb24 -y raw_video.raw'
  
  2. Converting a raw video file to mpeg 4:
    
    'avconv -f rawvideo -pix_fmt rgb24 -video_size 640x360 -framerate
          30 -i output_video.raw -pix_fmt yuv420p output_video.mp4'
          
The program works with RGB images of 3 bytes (rgb24). It is necessary to
know the dimensions of the original video and its framerate. The 
accompanying script 'cmdline_execute.sh' facilitates the task of converting
the videos to/from raw data and calling the midway program.

Usage instructions of the estadeo program:

  Usage: estadeo raw_input_video width height nframes [OPTIONS]

  Video stabilization program:
  'raw_input_video' is a video file in raw format (rgb24).
  'width' is the width of the images of the video in pixels.
  'height' is the height of the images of the video in pixels.
  'nframes' is the number of frames in the video.

  OPTIONS:
  --------
   -o name  output video name to write the computed raw video
              default value 'output_video.raw'
   -m N     motion estimation method
              0-direct method; 1-feature based method
              default value 1
   -t N     transformation type to be computed:
              2-traslation; 3-Euclidean transform;
              4-similarity; 6-affinity; 8-homography
              default value 8
   -s N     motion smoothing strategy:
              0-pure composition;
              1-compositional smoothing;
              2-compositional local smoothing; 
              3-local matrix-based smoothing;
              4-local point-based smoothing
              default value 3
   -r N     smoothing radius
              default value 30
   -b N     type of boundary condition: 
              0-constant; 1-neumann; 2-dirichlet
              default value 2
   -p N     video postprocessing 
              0-no postprocessing; 1-crop&zoom
              default value 0
   -w name  write transformations to file
   -l name  load transformations from file
   -v       switch on verbose mode 

   
Usage examples:

  1.Directly using the estadeo program:
    
  > avconv -i data/walk.mp4 -f rawvideo -pix_fmt rgb24 -y data/video.raw 
  > bin/estadeo data/video.raw 350 622 203 -v -o data/outvideo.raw -m 1 -s 3 \
           -r 50 -w data/transform.mat
  > avconv -f rawvideo -pix_fmt rgb24 -video_size 350x622 -framerate 30/1 -i \
           data/outvideo.raw -pix_fmt yuv420p -y data/stabilized.mp4
  
  2.Using the script:
    
   > bin/cmdline_execute.sh data/walk.mp4 data/outvideo.mp4 30 100 1 2 0.0000001 transforms.mat '-m 3 -b 1 -p 0 -online -t 4'
   

*************
LIST OF FILES
*************
main.cpp: Main function to be called from the command line. It reads and check 
  the parameters and call the estadeo program
estadeo.cpp: Implements the estadeo video stabilization
motion_smoothing.cpp: Implements the transformation smoothing strategies
video_cooling.cpp: Methods to improve the video after stabilization
utils.cpp: Functions to convert rgb videos to grayscale and add noise 
cmline_execute.sh: Script to be executed from the command line that facilitates
  the process of converting videos to/from raw data and calling the estadeo 
  algorithm

Complementary programs:
direct_method: a directory with the code of parametric motion estimation. 
  This is an implementation of the inverse compositional algorithm that is 
  available at http://www.ipol.im/pub/pre/153/
feature_based_method: a directory for feature-based image matching. This is
  the code from the ORSA method, http://www.ipol.im/pub/art/2012/mmm-oh/
generate_output.cpp: Program to create videos of the histograms of the input
  and output videos (used for the online demo only)
estadeo.sh: Script used from the IPOL demo to facilitate the process of 
  converting videos to/from raw data and calling the estadeo method
