#ifndef VBM3D_H_INCLUDED
#define VBM3D_H_INCLUDED

#include <fftw3.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"

//#define OPTICALFLOW

/** ------------------ **/
/** - Main functions - **/
/** ------------------ **/
//! Main function
int run_vbm3d(
    const float sigma
,   std::vector<float*> &buffer_input
,   std::vector<float*> &buffer_basic
,   float* final_estimate
,   const int w
,   const int h
,   const int d
,   const Parameters& prms_1
,   const Parameters& prms_2
,   const int index
,   const int size_buffer
,   std::vector<float>& kaiser_window_1
,   std::vector<float>& coef_norm_1
,   std::vector<float>& coef_norm_inv_1
,   std::vector<float>& kaiser_window_2
,   std::vector<float>& coef_norm_2
,   std::vector<float>& coef_norm_inv_2
,   std::vector<float>& lpd
,   std::vector<float>& hpd
,   std::vector<float>& lpr
,   std::vector<float>& hpr
);

//! 1st step of VBM3D
void vbm3d_1st_step(
    const float sigma
,   std::vector<float*>& buffer
,   float* basic_estimate
,   float* denominator
,   const int w
,   const int h
,   const int d
,   const Parameters& prms
,   fftwf_plan *  plan_2d
,   fftwf_plan *  plan_2d_inv
,   const int index
,   const int size_buffer
,   std::vector<float>& kaiser_window
,   std::vector<float>& coef_norm
,   std::vector<float>& coef_norm_inv
,   std::vector<float>& lpd
,   std::vector<float>& hpd
,   std::vector<float>& lpr
,   std::vector<float>& hpr
);

//! 2nd step of VBM3D
void vbm3d_2nd_step(
    const float sigma
,   std::vector<float*>& buffer_input
,   std::vector<float*>& buffer_basic
,   float* final_estimate
,   float* denominator
,   const int w
,   const int h
,   const int d
,   const Parameters& prms
,   fftwf_plan *  plan_2d
,   fftwf_plan *  plan_2d_inv
,   const int index
,   const int size_buffer
,   std::vector<float>& kaiser_window
,   std::vector<float>& coef_norm
,   std::vector<float>& coef_norm_inv
,   std::vector<float>& lpd
,   std::vector<float>& hpd
,   std::vector<float>& lpr
,   std::vector<float>& hpr
);

//! Process 2D dct of a group of patches
void dct_2d_process(
    std::vector<float> &DCT_table_2D
,   std::vector<float*> const& buffer
,   const int w
,   const int d
,   std::vector<std::pair<unsigned,unsigned> > const& patch_table 
,   fftwf_plan * plan
,   const unsigned kHW
,   std::vector<float> const& coef_norm
);

int computeSimilarPatches(
	std::vector<float>& output
,	std::vector<std::pair<unsigned,unsigned> >& index
,   int idx_curr_frame
,	unsigned pidx
, 	std::vector<float*>& frames
,   const int w
,   const int h
,   const int c
,   const int size_buffer
,	const Parameters& prms
);

//! Process 2D bior1.5 transform of a group of patches
void bior_2d_process(
    std::vector<float> &bior_table_2D
,   std::vector<float*> const& buffer
,   const int w
,   const int d
,   std::vector<std::pair<unsigned,unsigned> > const& patch_table
,   const unsigned kHW
,   std::vector<float> &lpd
,   std::vector<float> &hpd
);

void dct_2d_inv(
    std::vector<float> &group_3D_table
,   const unsigned kHW
,   const unsigned N
,   std::vector<float> const& coef_norm_inv
,   fftwf_plan * plan
);

void bior_2d_inv(
    std::vector<float> &group_3D_table
,   const unsigned kHW
,   std::vector<float> const& lpr
,   std::vector<float> const& hpr
);

//! HT filtering using Welsh-Hadamard transform (do only
//! third dimension transform, Hard Thresholding
//! and inverse Hadamard transform)
void ht_filtering_hadamard(
    std::vector<float> &group_3D
,   std::vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kHard
,   const unsigned chnls
,   float sigma
,   const float lambdaThr3D
,   std::vector<float> &weight_table
);

//! HT filtering using Haar transform (do only
//! third dimension transform, Hard Thresholding
//! and inverse Hadamard transform)
void ht_filtering_haar(
    std::vector<float> &group_3D
,   std::vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kHard
,   const unsigned chnls
,   float sigma
,   const float lambdaThr3D
,   std::vector<float> &weight_table
);

//! Wiener filtering using Welsh-Hadamard transform
void wiener_filtering_hadamard(
    std::vector<float> &group_3D_img
,   std::vector<float> &group_3D_est
,   std::vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kWien
,   const unsigned chnls
,   float sigma
,   std::vector<float> &weight_table
);

//! Wiener filtering using Haar transform
void wiener_filtering_haar(
    std::vector<float> &group_3D_img
,   std::vector<float> &group_3D_est
,   std::vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kWien
,   const unsigned chnls
,   float sigma
,   std::vector<float> &weight_table
);

/** ---------------------------------- **/
/** - Preprocessing / Postprocessing - **/
/** ---------------------------------- **/
//! Preprocess coefficients of the Kaiser window and normalization coef for the DCT
void preProcess(
    std::vector<float> &kaiserWindow
,   std::vector<float> &coef_norm
,   std::vector<float> &coef_norm_inv
,   const unsigned kHW
);

#endif // VBM3D_H_INCLUDED
