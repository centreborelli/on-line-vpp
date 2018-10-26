// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2017, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#include <stdio.h>


/**
  *
  *  Function to read a video in raw data
  * 
**/
size_t read_video(
  char *name,       //file name
  unsigned char *I, //video to read
  int size          //size of the video
) 
{
  FILE* file=fopen(name, "rb");
  if(file != NULL)
  {
    size_t result=fread(I, sizeof(unsigned char), size, file);
    fclose(file);
    return result;
  }
  else return 0;
}


/**
  *
  *  Function to read a video in raw data
  * 
**/
size_t write_video(
  char *name,       //file name
  unsigned char *I, //video to write
  int size          //size of the video
) 
{
  FILE* file=fopen(name, "wb");
  if(file != NULL)
  {
    size_t result=fwrite(I, sizeof(unsigned char), size, file);
    fclose(file);
    return result;
  }
  else return 0;
}


/**
  *
  *  Function to save transformations to a file
  * 
**/
void save_transform(
  char  *name,   //file name
  float *H,      //transformation
  int   nparams  //number of parameters of the transformations
)
{
  FILE *fd=fopen(name,"a");
  if(fd!=NULL)
  {
    for(int j=0;j<nparams;j++) fprintf(fd, "%.10f ", H[j]);
    fprintf(fd, "\n");
    fclose(fd);
  }
}

