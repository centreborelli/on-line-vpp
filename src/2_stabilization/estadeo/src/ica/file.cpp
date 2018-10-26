// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#include <stdio.h>

#include "file.h"


/**
 *
 *  Function to read the parameters in ascii format
 *  It reads a header with: nparams nx ny
 *  Then it reads the parameters for each pixel
 *
 */
void read
(
  char *file,   //input file name
  float **p,    //parameters to be read
  int &nparams, //number of parameters
  int &nx,      //number of columns 
  int &ny       //number of rows 
)
{
  FILE *fd=fopen(file,"r");
  if(fd!=NULL)
  {
    int result=fscanf(fd,"%d %d %d",&nparams, &nx, &ny);
    if(result==3)
    {
      *p=new float[nx*ny*nparams];
      for(int i=0; i<nx*ny; i++)
      {
        for(int j=0; j<nparams; j++)
          result=fscanf(fd,"%f", &((*p)[i*nparams+j]));
      }
    }
    fclose(fd);
  }
}


/**
 *
 *  Function to save the parameters in ascii format
 *  It creates a header with: nparams nx ny
 *  Then it stores the parameters for each pixel
 *
 */
void save
(
  char *file,  //output file name 
  float *p,   //parameters to be saved
  int nparams, //number of parameters
  int nx,      //number of columns 
  int ny       //number of rows 
) 
{
  
  FILE *fd=fopen(file,"w");
  fprintf(fd,"%d %d %d\n", nparams, nx, ny);
  for(int i=0; i<nx*ny; i++)
  {
    for(int j=0; j<nparams; j++)
      fprintf(fd,"%f ", p[i*nparams+j]);
    fprintf(fd,"\n");
  }
  fclose(fd);
}


/**
 *
 *  Function to read the parameters in ascii format
 *  It reads a header with: nparams nx ny
 *  Then it reads the parameters 
 *
 */
void read
(
  const char *file, //input file name
  float **p,       //parameters to be read
  int &nparams      //number of parameters
)
{
  FILE *fd=fopen(file,"r");
  *p=NULL;
  if(fd!=NULL)
  {
    int result=fscanf(fd,"%d",&nparams);
    if(result==1)
    {
      *p=new float[nparams];
      for(int j=0; j<nparams; j++)
        result=fscanf(fd,"%f", &((*p)[j]));
    }
    fclose(fd);
  }
}

/**
 *
 *  Function to save the parameters in ascii format
 *  It creates a header with: nparams nx ny
 *  Then it stores the parameters 
 *
 */
void save
(
  const char *file, //output file name 
  float *p,        //parameters to be saved
  int nparams       //number of parameters
) 
{
  FILE *fd=fopen(file,"w");
  fprintf(fd,"%d\n", nparams);
  for(int j=0; j<nparams; j++)
    fprintf(fd,"%f ", p[j]);
  fprintf(fd,"\n");
  fclose(fd);
}
