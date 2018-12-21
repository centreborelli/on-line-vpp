/*
 * Copyright (c) 2018, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>
#include <algorithm>

#include "vbm3d.h"
#include "vpp.h"
#include "Utilities.h"
#include "cmd_option.h"
#include "lib_transforms.h"

#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3
#define DCT       4
#define BIOR      5
#define HADAMARD  6
#define HAAR      7
#define NONE      8

using namespace std;

/**
 * @brief Initialize parameters whose default value is constant.
 *
 **/
void initConstantParams_1(
	Parameters& prms
,	const int k
,	const int Rf
,	const int Rr
,	const int Rp
,	const int Nl
,	const int st
,	const int Nmax
,	const int d
,	const float lambda
,	const unsigned T_2D
,	const unsigned T_3D
){
	const float n = 1./255./255.;
	prms.k    = (k    < 0) ? 8 : k;
	prms.Rf   = (Rf   < 0) ? 4 : Rf;
	prms.Rr   = (Rr   < 0) ? 7 : Rr;
	prms.Rp   = (Rp   < 0) ? 5 : Rp;
	prms.Nl   = (Nl   < 0) ? 2 : Nl;
	prms.d    = (d    < 0) ? 7.*7.*n : d*d*n;
	prms.st   = (st   < 0) ? 6 : st;
	prms.Nmax = (Nmax < 0) ? 8 : Nmax;
	prms.lambda = (lambda < 0) ? 2.7f : lambda;
	prms.T_2D = (T_2D == NONE) ? BIOR : T_2D;
	prms.T_3D = (T_3D == NONE) ? HAAR : T_3D;
}

void initConstantParams_2(
	Parameters& prms
,	const int Rf
,	const int Rr
,	const int Rp
,	const int Nl
,	const int st
,	const int Nmax
,	const int d
,	const unsigned T_2D
,	const unsigned T_3D
){
	const float n = 1./255./255.;
	prms.Rf   = (Rf   < 0) ? 4 : Rf;
	prms.Rr   = (Rr   < 0) ? 7 : Rr;
	prms.Rp   = (Rp   < 0) ? 5 : Rp;
	prms.Nl   = (Nl   < 0) ? 2 : Nl;
	prms.d    = (d    < 0) ? 3.*3.*n : d*d*n;
	prms.st   = (st   < 0) ? 4 : st;
	prms.Nmax = (Nmax < 0) ? 8 : Nmax;
	prms.T_2D = (T_2D == NONE) ? DCT  : T_2D;
	prms.T_3D = (T_3D == NONE) ? HAAR : T_3D;
}

/**
 * @brief Initialize parameters whose default value vary with sigma
 *
 **/
void initSigmaParams_1(
	Parameters& prms
,	const float tau
,	const float sigma
){
	const float n = 1./255./255.;
	prms.tau = (tau < 0) ? ( (sigma > 30) ? 4500 : 3000 )*n : tau;
}

void initSigmaParams_2(
	Parameters& prms
,	const int k
,	const float tau
,	const float sigma
){
	const float n = 1./255./255.;
	prms.k   = (k   < 0) ? ( (sigma > 30) ? 8 : 7 ) : k;
	prms.tau = (tau < 0) ? ( (sigma > 30) ? 3000 : 1500 )*n : tau;
}


/**
 * @file   main.cpp
 * @brief  Main executable file. Do not use lib_fftw to
 *         process DCT.
 */
int main(int argc, char **argv)
{
	clo_usage("Causal version of the VBM3D video denoising method");
	clo_help(" NOTE: Input (<) and output (>) sequences have to be vpp pipes.\n");

	using std::string;
	const string  input_path = clo_option("-i", "-", "< input pipe");
	const string  final_path = clo_option("-o", "-", "> output pipe");
	const string  sigma_path = clo_option("-s", "", "< noise sigma pipe");

	//! VBM3D parameters
	const int   Rf      = clo_option("-Rf"     ,-1 , "< number of previous frames");
	const int   Rr      = clo_option("-Rr"     ,-1 , "< size of searach region in current frame");
	const int   Rp      = clo_option("-Rp"     ,-1 , "< size of local search regions in previous frames");
	const int   Nl      = clo_option("-Nl"     ,-1 , "< number of candidates per local search region");
	const int   kH      = clo_option("-kH"     ,-1 , "< patch size");
	const int   kW      = clo_option("-kW"     ,-1 , "< patch size");
	const int   stH     = clo_option("-stH"    ,-1 , "< step to skip patches");
	const int   stW     = clo_option("-stW"    ,-1 , "< step to skip patches");
	const int   NmaxH   = clo_option("-NmaxH"  ,-1 , "< maximum number of similar patches");
	const int   NmaxW   = clo_option("-NmaxW"  ,-1 , "< maximum number of similar patches");
	const int   dH      = clo_option("-dH"     ,-1 , "< d offset in distance");
	const int   dW      = clo_option("-dW"     ,-1 , "< d offset in distance");
	const float tauH    = clo_option("-tauH"   ,-1., "< distance threshold");
	const float tauW    = clo_option("-tauW"   ,-1., "< distance threshold");
	const float lambdaH = clo_option("-lambdaH",-1., "< parameter for the thresholding operator");
	const unsigned color_space  =  (unsigned) clo_option("-color", 0 , "< color space");
	const unsigned T2DH  = (unsigned) clo_option("-T2DH", NONE , "< Spatial transform for hard thresholding step");
	const unsigned T2DW  = (unsigned) clo_option("-T2DW", NONE , "< Spatial transform for Wiener filtering step");
	const unsigned T3DH  = (unsigned) clo_option("-T3DH", NONE , "< 3rd dimension transform for hard thresholding step");
	const unsigned T3DW  = (unsigned) clo_option("-T3DW", NONE , "< 3rd dimension transform for Wiener filtering step");

	//! Check inputs
	if (input_path == "")
		return fprintf(stderr, "%s: no input images.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]), 1;

	if (T2DH != NONE && T2DH != DCT && T2DH != BIOR)
		return fprintf(stderr, "Unknown T2D_H: %d for DCT %d for bi-orthogonal\n", DCT, BIOR), 1;

	if (T2DW != NONE && T2DW != DCT && T2DW != BIOR)
		return fprintf(stderr, "Unknown T2D_W: %d for DCT %d for bi-orthogonal\n", DCT, BIOR), 1;

	if (T3DH != NONE && T3DH != HAAR && T3DH != HADAMARD)
		return fprintf(stderr, "Unknown T3D_H: %d for HAAR %d for HADAMARD\n", HAAR, HADAMARD), 1;

	if (T3DW != NONE && T3DW != HAAR && T3DW != HADAMARD)
		return fprintf(stderr, "Unknown T3D_W: %d for HAAR %d for HADAMARD\n", HAAR, HADAMARD), 1;

	//! Initialize parameters independent of the noise level
	Parameters prms_1;
	Parameters prms_2;
	initConstantParams_1(prms_1, kH, Rf, Rr, Rp, Nl, stH, NmaxH, dH, lambdaH, T2DH, T3DH);
	initConstantParams_2(prms_2,     Rf, Rf, Rp, Nl, stW, NmaxW, dW,          T2DW, T3DW);

	//! Init Kaiser Window for first step (since it doesn't depend on sigma)
	int kH2 = prms_1.k*prms_1.k;
	vector<float> kaiser_window_1(kH2);
	vector<float> coef_norm_1(kH2);
	vector<float> coef_norm_inv_1(kH2);
	kaiserWindow(kaiser_window_1, coef_norm_1, coef_norm_inv_1, prms_1.k);

	//! Preprocessing of Bior table
	vector<float> lpd, hpd, lpr, hpr;
	bior15_coef(lpd, hpd, lpr, hpr);

	//! Init pipes
	int w,h,d;
	FILE* in = vpp_init_input(input_path.c_str(), &w, &h, &d);
	if (!in)
		return fprintf(stderr, "vbm3d: cannot initialize input '%s'\n", input_path.c_str()), 1;

	FILE* out = vpp_init_output(final_path.c_str(), w, h, d);
	if (!out)
		return fprintf(stderr, "vbm3d: cannot initialize output '%s'\n", final_path.c_str()), 1;

	int sw, sh, sd;
	FILE* sigma_in = vpp_init_input(sigma_path.c_str(), &sw, &sh, &sd);
	float sigma;
	if (!sigma_in)
		return fprintf(stderr, "vbm3d: cannot initialize sigma '%s'\n", sigma_path.c_str()), 1;

	if (sw != 2 || sh != 1 || sd != 1)
		return fprintf(stderr, "vbm3d: invalid sigma stream size: %dx%dx%d\n", sw, sh, sd), 1;

	//! Initialize the buffers to store the previous frames
	vector<float*> buffer_input(prms_1.Rf, NULL);
	vector<float*> buffer_basic(prms_2.Rf, NULL);
	for(int i = 0; i < prms_1.Rf; ++i)
	{
		buffer_input[i] = (float*) malloc(w*h*d*sizeof(*(buffer_input[i])));
		buffer_basic[i] = (float*) malloc(w*h*d*sizeof(*(buffer_basic[i])));
	}
	//! Buffer for denoised frame
	float* final_estimate;
	final_estimate = (float*) malloc(w*h*d*sizeof*final_estimate);

	//! Process the frames
	int size_buffer = 0;
	int index = 0;
	while (vpp_read_frame(in, buffer_input[index], w, h, d))
	{
		// Update buffer size
		size_buffer = std::min(size_buffer+1, (int)prms_1.Rf);

		// Read noise level
		float sigmadata[2];
		vpp_read_frame(sigma_in, sigmadata, sw, sh, sd);
		sigma = sigmadata[1];

		//! Initialize parameters independent of the noise level
		initSigmaParams_1(prms_1,     tauH, sigma);
		initSigmaParams_2(prms_2, kW, tauW, sigma);

		//! Init Kaiser Window for second step
		int kW2 = prms_2.k*prms_2.k;
		vector<float> kaiser_window_2(kW2);
		vector<float> coef_norm_2(kW2);
		vector<float> coef_norm_inv_2(kW2);
		kaiserWindow(kaiser_window_2, coef_norm_2, coef_norm_inv_2, prms_2.k);

		// Change colorspace (RGB to OPP)
		if(color_space == 0)
			transformColorSpace(buffer_input[index], w, h, d, true);

		run_vbm3d(sigma, buffer_input, buffer_basic, final_estimate, w, h, d,
				prms_1, prms_2, index, size_buffer, kaiser_window_1, coef_norm_1,
				coef_norm_inv_1, kaiser_window_2, coef_norm_2, coef_norm_inv_2,
				lpd, hpd, lpr, hpr);

		// Inverse the colorspace for the output (OPP to RGB)
		if(color_space == 0)
			transformColorSpace(final_estimate, w, h, d, false);

		// send the frame to the next step
		if (!vpp_write_frame(out, final_estimate, w, h, d))
			break;

		index = (index+1) % prms_1.Rf;
	}

	for(int i = 0; i < prms_1.Rf; ++i)
	{
		free(buffer_input[i]);
		free(buffer_basic[i]);
	}
	free(final_estimate);
	fclose(in);
	fclose(out);

	return EXIT_SUCCESS;
}
