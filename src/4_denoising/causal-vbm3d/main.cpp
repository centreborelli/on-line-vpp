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
,	const int Nf
,	const int Ns
,	const int Npr
,	const int Nb
,	const int p
,	const int N
,	const int d
,	const float lambda3D
,	const unsigned T_2D
,	const unsigned T_3D
){
	const float n = 1./255./255.;
	prms.k   = (k   < 0) ? 8 : k;
	prms.Nf  = (Nf  < 0) ? 4 : Nf;
	prms.Ns  = (Ns  < 0) ? 7 : Ns;
	prms.Npr = (Npr < 0) ? 5 : Npr;
	prms.Nb  = (Nb  < 0) ? 2 : Nb;
	prms.d   = (d   < 0) ? 7.*7.*n : d*d*n;
	prms.p   = (p   < 0) ? 6 : p;
	prms.N   = (N   < 0) ? 8 : N;
	prms.lambda3D = (lambda3D < 0) ? 2.7f : lambda3D;
	prms.T_2D = (T_2D == NONE) ? BIOR : T_2D;
	prms.T_3D = (T_3D == NONE) ? HAAR : T_3D;
}

void initConstantParams_2(
	Parameters& prms
,	const int Nf
,	const int Ns
,	const int Npr
,	const int Nb
,	const int p
,	const int N
,	const int d
,	const unsigned T_2D
,	const unsigned T_3D
){
	const float n = 1./255./255.;
	prms.Nf  = (Nf  < 0) ? 4 : Nf;
	prms.Ns  = (Ns  < 0) ? 7 : Ns;
	prms.Npr = (Npr < 0) ? 5 : Npr;
	prms.Nb  = (Nb  < 0) ? 2 : Nb;
	prms.d   = (d   < 0) ? 3.*3.*n : d*d*n;
	prms.p   = (p   < 0) ? 4 : p;
	prms.N   = (N   < 0) ? 8 : N;
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
	const int   kHard    = clo_option("-kH"   , -1 , "< patch size");
	const int   NfHard   = clo_option("-NfH"  , -1 , "< number of previous frames");
	const int   NsHard   = clo_option("-NsH"  , -1 , "< size of region if current frame");
	const int   NprHard  = clo_option("-NprH" , -1 , "< size of local search regions in previous frames");
	const int   NbHard   = clo_option("-NbH"  , -1 , "< number of candidates per local search region");
	const int   pHard    = clo_option("-pH"   , -1 , "< step to skip patches");
	const int   NHard    = clo_option("-NH"   , -1 , "< maximum number of similar patches");
	const int   dHard    = clo_option("-dH"   , -1 , "< d offset in distance");
	const float tauHard  = clo_option("-tauH" , -1., "< distance threshold");
	const float lambda3D = clo_option("-lambda",-1., "< parameter for the thresholding operator");
	const int   kWien    = clo_option("-kW"   , -1 , "< patch size");
	const int   NfWien   = clo_option("-NfW"  , -1 , "< number of previous frames");
	const int   NsWien   = clo_option("-NsW"  , -1 , "< size of region if current frame");
	const int   NprWien  = clo_option("-NprW" , -1 , "< size of local search regions in previous frames");
	const int   NbWien   = clo_option("-NbW"  , -1 , "< number of candidates per local search region");
	const int   pWien    = clo_option("-pW"   , -1 , "< step to skip patches");
	const int   NWien    = clo_option("-NW"   , -1 , "< maximum number of similar patches");
	const int   dWien    = clo_option("-dW"   , -1 , "< d offset in distance");
	const float tauWien  = clo_option("-tauW" , -1., "< distance threshold");
	const unsigned color_space  =  (unsigned) clo_option("-color", 0 , "< color space");
	const unsigned T_2D_hard  = (unsigned) clo_option("-T2dh", NONE , "< Spatial transform for hard thresholding step");
	const unsigned T_2D_wien  = (unsigned) clo_option("-T2dw", NONE , "< Spatial transform for Wiener filtering step");
	const unsigned T_3D_hard  = (unsigned) clo_option("-T3dh", NONE , "< 3rd dimension transform for hard thresholding step");
	const unsigned T_3D_wien  = (unsigned) clo_option("-T3dw", NONE , "< 3rd dimension transform for Wiener filtering step");

	//! Check inputs
	if (input_path == "")
		return fprintf(stderr, "%s: no input images.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]), 1;

	if (T_2D_hard != NONE && T_2D_hard != DCT && T_2D_hard != BIOR)
		return fprintf(stderr, "Unknown T2D_H: %d for DCT %d for bi-orthogonal\n", DCT, BIOR), 1;

	if (T_2D_wien != NONE && T_2D_wien != DCT && T_2D_wien != BIOR)
		return fprintf(stderr, "Unknown T2D_W: %d for DCT %d for bi-orthogonal\n", DCT, BIOR), 1;

	if (T_3D_hard != NONE && T_3D_hard != HAAR && T_3D_hard != HADAMARD)
		return fprintf(stderr, "Unknown T3D_H: %d for HAAR %d for HADAMARD\n", HAAR, HADAMARD), 1;

	if (T_3D_wien != NONE && T_3D_wien != HAAR && T_3D_wien != HADAMARD)
		return fprintf(stderr, "Unknown T3D_W: %d for HAAR %d for HADAMARD\n", HAAR, HADAMARD), 1;

	//! Initialize parameters independent of the noise level
	Parameters prms_1;
	Parameters prms_2;
	initConstantParams_1(prms_1, kHard, NfHard, NsHard, NprHard, NbHard, pHard, NHard, dHard, lambda3D, T_2D_hard, T_3D_hard);
	initConstantParams_2(prms_2,        NfWien, NsWien, NprWien, NbWien, pWien, NWien, dWien,           T_2D_wien, T_3D_wien);

	// Force both buffer to have the same size for now
	prms_2.Nf = prms_1.Nf;

	//! Init Kaiser Window for first step (since it doesn't depend on sigma)
	int kHard_2 = prms_1.k*prms_1.k;
	vector<float> kaiser_window_1(kHard_2);
	vector<float> coef_norm_1(kHard_2);
	vector<float> coef_norm_inv_1(kHard_2);
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
	vector<float*> buffer_input(prms_1.Nf, NULL);
	vector<float*> buffer_basic(prms_2.Nf, NULL);
	for(int i = 0; i < prms_1.Nf; ++i)
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
		size_buffer = std::min(size_buffer+1, (int)prms_1.Nf);

		// Read noise level
		float sigmadata[2];
		vpp_read_frame(sigma_in, sigmadata, sw, sh, sd);
		sigma = sigmadata[1];

		//! Initialize parameters independent of the noise level
		initSigmaParams_1(prms_1,        tauHard, sigma);
		initSigmaParams_2(prms_2, kWien, tauWien, sigma);

		//! Init Kaiser Window for second step
		int kWien_2 = prms_2.k*prms_2.k;
		vector<float> kaiser_window_2(kWien_2);
		vector<float> coef_norm_2(kWien_2);
		vector<float> coef_norm_inv_2(kWien_2);
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

		index = (index+1) % prms_1.Nf;
	}

	for(int i = 0; i < prms_1.Nf; ++i)
	{
		free(buffer_input[i]);
		free(buffer_basic[i]);
	}
	free(final_estimate);
	fclose(in);
	fclose(out);

	return EXIT_SUCCESS;
}
