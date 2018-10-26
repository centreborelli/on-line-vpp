// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2015, Javier Sánchez Pérez <jsanchez@ulpgc.es>
// All rights reserved.


#ifndef FILE_H
#define FILE_H


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
  float **p,   //parameters to be read
  int &nparams, //number of parameters
  int &nx,      //number of columns 
  int &ny       //number of rows 
);

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
);

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
);


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
);

#endif
