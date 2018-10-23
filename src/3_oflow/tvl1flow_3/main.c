
// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2011, Javier Sánchez Pérez <jsanchez@dis.ulpgc.es>
// All rights reserved.

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#ifndef DISABLE_OMP
#include <omp.h>
#endif//DISABLE_OMP

#include "vpp.h"

#include "tvl1flow_lib.h"


#define PAR_DEFAULT_NPROC   0
#define PAR_DEFAULT_TAU     0.25
#define PAR_DEFAULT_LAMBDA  0.15
#define PAR_DEFAULT_THETA   0.3
#define PAR_DEFAULT_NSCALES 100
#define PAR_DEFAULT_FSCALE  0
#define PAR_DEFAULT_ZFACTOR 0.5
#define PAR_DEFAULT_NWARPS  5
#define PAR_DEFAULT_EPSILON 0.01
#define PAR_DEFAULT_VERBOSE 0


/**
 *
 *  Main program:
 *   This program reads the following parameters from the console and
 *   then computes the optical flow:
 *   -nprocs      number of threads to use (OpenMP library)
 *   -vid_path    pipe in
 *   -out_path    pipe for the output flow field
 *   -tau         time step in the numerical scheme
 *   -lambda      data term weight parameter
 *   -theta       tightness parameter
 *   -nscales     number of scales in the pyramidal structure
 *   -zfactor     downsampling factor for creating the scales
 *   -nwarps      number of warps per scales
 *   -epsilon     stopping criterion threshold for the iterative process
 *   -verbose     switch on/off messages
 *
 */
int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: %s vid  out ["
		//                       0 1    2
		"nproc tau lambda theta nscales fscale zfactor nwarps epsilon "
		//  3  4   5      6     7       8       9      10     11 
		"verbose]\n", *argv);
		// 12
		return EXIT_FAILURE;
	}

	//read the parameters
	int i = 1;
	char* vid_path  = argv[i]; i++;
	char* flow_path = argv[i]; i++;
	int   nproc   = (argc>i)? atoi(argv[i]): PAR_DEFAULT_NPROC;   i++;
	float tau     = (argc>i)? atof(argv[i]): PAR_DEFAULT_TAU;     i++;
	float lambda  = (argc>i)? atof(argv[i]): PAR_DEFAULT_LAMBDA;  i++;
	float theta   = (argc>i)? atof(argv[i]): PAR_DEFAULT_THETA;   i++;
	int   nscales = (argc>i)? atoi(argv[i]): PAR_DEFAULT_NSCALES; i++;
	int   fscale  = (argc>i)? atoi(argv[i]): PAR_DEFAULT_FSCALE ; i++;
	float zfactor = (argc>i)? atof(argv[i]): PAR_DEFAULT_ZFACTOR; i++;
	int   nwarps  = (argc>i)? atoi(argv[i]): PAR_DEFAULT_NWARPS;  i++;
	float epsilon = (argc>i)? atof(argv[i]): PAR_DEFAULT_EPSILON; i++;
	int   verbose = (argc>i)? atoi(argv[i]): PAR_DEFAULT_VERBOSE; i++;

	//check parameters
	if (nproc < 0) {
		nproc = PAR_DEFAULT_NPROC;
		if (verbose) fprintf(stderr, "warning: "
				"nproc changed to %d\n", nproc);
	}
	if (tau <= 0 || tau > 0.25) {
		tau = PAR_DEFAULT_TAU;
		if (verbose) fprintf(stderr, "warning: "
				"tau changed to %g\n", tau);
	}
	if (lambda <= 0) {
		lambda = PAR_DEFAULT_LAMBDA;
		if (verbose) fprintf(stderr, "warning: "
				"lambda changed to %g\n", lambda);
	}
	if (theta <= 0) {
		theta = PAR_DEFAULT_THETA;
		if (verbose) fprintf(stderr, "warning: "
				"theta changed to %g\n", theta);
	}
	if (nscales <= 0) {
		nscales = PAR_DEFAULT_NSCALES;
		if (verbose) fprintf(stderr, "warning: "
				"nscales changed to %d\n", nscales);
	}
	if (zfactor <= 0 || zfactor >= 1) {
		zfactor = PAR_DEFAULT_ZFACTOR;
		if (verbose) fprintf(stderr, "warning: "
				"zfactor changed to %g\n", zfactor);
	}
	if (nwarps <= 0) {
		nwarps = PAR_DEFAULT_NWARPS;
		if (verbose) fprintf(stderr, "warning: "
				"nwarps changed to %d\n", nwarps);
	}
	if (epsilon <= 0) {
		epsilon = PAR_DEFAULT_EPSILON;
		if (verbose) fprintf(stderr, "warning: "
				"epsilon changed to %f\n", epsilon);
	}

#ifndef DISABLE_OMP
	if (nproc > 0)
		omp_set_num_threads(nproc);
#endif//DISABLE_OMP

	// read the input images
	int nx, ny, c;

	FILE* vid_in = vpp_init_input(vid_path, &nx, &ny, &c);
	if (!vid_in)
	{
		fprintf(stderr, "%s: cannot initialize input '%s'\n", argv[0], vid_path);
		return EXIT_FAILURE;
	}

	FILE* flow_out = vpp_init_output(flow_path, nx, ny, 2);
	if (!flow_out)
	{
		fprintf(stderr, "%s: cannot initialize output '%s'\n", argv[0], flow_path);
		return EXIT_FAILURE;
	}

	// allocate memory for frames and optical flow
	float *I0 = malloc(nx*ny*c*sizeof(float));
	float *I1 = malloc(nx*ny*c*sizeof(float));
	float *u = malloc(2*nx*ny*sizeof(float));
	float *v = u + nx*ny;;

	// FIXME: for conversion between split layout to channels layout
	float *uv = malloc(2*nx*ny*sizeof(float));

	// read first frame
	vpp_read_frame(vid_in, I0, nx, ny, c);

	//read the images and compute the optical flow
//	int f = 1;
	while (vpp_read_frame(vid_in, I1, nx, ny, c))
	{
		//Set the number of scales according to the size of the
		//images.  The value N is computed to assure that the smaller
		//images of the pyramid don't have a size smaller than 16x16
		const float N = 1 + log(hypot(nx, ny)/16.0) / log(1/zfactor);
		if (N < nscales) nscales = N;
		if (nscales < fscale) fscale = nscales;

		if (verbose)
			fprintf(stderr,
				"nproc=%d tau=%f lambda=%f theta=%f nscales=%d "
				"zfactor=%f nwarps=%d epsilon=%g\n",
				nproc, tau, lambda, theta, nscales,
				zfactor, nwarps, epsilon);

		//compute the optical flow
#define BACKWARDS_FLOW
#ifdef BACKWARDS_FLOW
		Dual_TVL1_optic_flow_multiscale(
				I1, I0, u, v, nx, ny, tau, lambda, theta,
				nscales, fscale, zfactor, nwarps, epsilon, verbose
		);
#else
		Dual_TVL1_optic_flow_multiscale(
				I0, I1, u, v, nx, ny, tau, lambda, theta,
				nscales, fscale, zfactor, nwarps, epsilon, verbose
		);
#endif

		//convert to vector layout
		for (int i = 0; i < nx*ny; ++i)
		{
			uv[2*i + 0] = u[i];
			uv[2*i + 1] = v[i];
		}

		//save the optical flow
		vpp_write_frame(flow_out, uv, nx, ny, 2);

		//swap buffers
		float *tmp = I0; I0 = I1; I1 = tmp;
//		f++;
	}

	//delete allocated memory
	free(u);
	free(uv);
	free(I0);
	free(I1);

	return EXIT_SUCCESS;
}
