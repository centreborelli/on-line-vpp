/*
 * Modified work copyright (c) 2018, Thibaud Ehret <ehret.thibaud@gmail.com>
 * Original work copyright (c) 2011, Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @file vbm3d.cpp
 * @brief VBM3D denoising functions
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include <iostream>
#include <algorithm>
#include <math.h>
#include <unordered_map>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "vbm3d.h"
#include "Utilities.h"
#include "lib_transforms.h"

#define SQRT2     1.414213562373095
#define SQRT2_INV 0.7071067811865475
#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3
#define DCT       4
#define BIOR      5
#define HADAMARD  6
#define HAAR      7

#ifdef _OPENMP
#include <omp.h>
#endif

/* 
 * In order to reproduce the original VBM3D the DC coefficients are
 * not thresholded (DCTHRESH commented) but are filtered using Wiener
 * (DCWIENER uncommented), MTRICK activates undocumented tricks from 
 * Marc Lebrun's implementation of BM3D available in IPOL 
 * http://www.ipol.im/pub/art/2012/l-bm3d/ 
 */

//#define DCTHRESH
#define DCWIENER
//#define MTRICK

using namespace std;

bool ComparaisonFirst(pair<float,unsigned> pair1, pair<float,unsigned> pair2)
{
	return pair1.first < pair2.first;
}

/** ----------------- **/
/** - Main function - **/
/** ----------------- **/
/**
 * @brief run VBM3D process.
 *
 * @param sigma: value of assumed noise of the noisy video;
 * @param buffer_input: vector containing the noisy frames;
 * @param buffer_basic: vector containing the basic estimations;
 * @param final_estimate: will be the denoised frame;
 * @param w: width of a frame;
 * @param h: height of a frame;
 * @param d: number of channels of a frame;
 * @param prms_1: parameters for the first step;
 * @param prms_2: parameters for the second step;
 * @param index: index of the frame in the buffer;
 * @param size_buffer: current size of the buffer;
 * @param kaiser_window_1: coefficient for the kaiser window for the first step;
 * @param coef_norm_1: normalization coeffs for the first step;
 * @param coef_norm_inv_1: inverse normalization coeffs for the first step;
 * @param kaiser_window_2: coefficient for the kaiser window for the second step;
 * @param coef_norm_2: normalization coeffs for the second step;
 * @param coef_norm_inv_2: inverse normalization coeffs for the second step;
 * @param lpd, hpd, lpr, hpr: bior basis;
 *
 * @return EXIT_FAILURE if color_space has not expected
 *         type, otherwise return EXIT_SUCCESS.
 **/
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
,   vector<float>& kaiser_window_1
,   vector<float>& coef_norm_1
,   vector<float>& coef_norm_inv_1
,   vector<float>& kaiser_window_2
,   vector<float>& coef_norm_2
,   vector<float>& coef_norm_inv_2
,   vector<float>& lpd
,   vector<float>& hpd
,   vector<float>& lpr
,   vector<float>& hpr
){
    float* denominator;
    denominator = (float*) malloc(w*h*d*sizeof*denominator);

	//! Allocate plan for FFTW library
	fftwf_plan plan_2d[1];
	fftwf_plan plan_2d_inv[1];

    //! Allocating Plan for FFTW process
    if (prms_1.T_2D == DCT)
    {
        const unsigned nb_cols = ind_size(0, w - prms_1.k, prms_1.p);
        allocate_plan_2d(&plan_2d[0], prms_1.k, FFTW_REDFT10,
                prms_1.N * d);
        allocate_plan_2d(&plan_2d_inv[0], prms_1.k, FFTW_REDFT01,
                prms_1.N * nb_cols * d);
    }

    //! Denoising, 1st Step
    vbm3d_1st_step(sigma, buffer_input, buffer_basic[index], denominator, w, h, d, prms_1,
            &plan_2d[0], &plan_2d_inv[0], index, size_buffer, kaiser_window_1, coef_norm_1,
            coef_norm_inv_1, lpd, hpd, lpr, hpr);

    //! Aggregation basic 1st Step
    for (unsigned k = 0; k < w*h*d; k++)
        buffer_basic[index][k] /= denominator[k];

    if(prms_2.k > 0)
    {
        //! Allocating Plan for FFTW process
        if (prms_2.T_2D == DCT)
        {
            const unsigned nb_cols = ind_size(0, (w - prms_2.k), prms_2.p);
            allocate_plan_2d(&plan_2d[0], prms_2.k, FFTW_REDFT10,
                    prms_2.N * d);
            allocate_plan_2d(&plan_2d_inv[0], prms_2.k, FFTW_REDFT01,
                    prms_2.N * nb_cols * d);
        }
        //! Denoising, 2nd Step
        vbm3d_2nd_step(sigma, buffer_input, buffer_basic, final_estimate, denominator, w, h, d, prms_2,
                &plan_2d[0], &plan_2d_inv[0], index, size_buffer, kaiser_window_2, coef_norm_2,
                coef_norm_inv_2, lpd, hpd, lpr, hpr);

        //! Aggregation 2nd Step
        for (unsigned k = 0; k < w*h*d; k++)
            final_estimate[k] /= denominator[k];
    }
    else
    {
        for (unsigned k = 0; k < w*h*d; k++)
            final_estimate[k] = buffer_basic[index][k];
    }

	return EXIT_SUCCESS;
}

/**
 * @brief Run the basic process of BM3D (1st step). The result
 *        is contained in basic_estimate.
 *
 * @param sigma: value of assumed noise of the video to denoise;
 * @param buffer: vector containing the noisy frames;
 * @param basic_estimate: will be the denoised frame;
 * @param denominator: will containes the aggregation coefficients;
 * @param w: width of a frame;
 * @param h: height of a frame;
 * @param d: number of channels of a frame;
 * @param prms: parameters;
 * @param plan_2d: fftw plan for the dct transform;
 * @param plan_2d_inv: fftw plan for the inverse dct transform;
 * @param index: index of the frame in the buffer;
 * @param size_buffer: current size of the buffer;
 * @param kaiser_window: coefficient for the kaiser window;
 * @param coef_norm: normalization coeffs;
 * @param coef_norm_inv: inverse normalization coeffs;
 * @param lpd, hpd, lpr, hpr: bior basis;
 *
 * @return none.
 **/
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
,   vector<float>& kaiser_window
,   vector<float>& coef_norm
,   vector<float>& coef_norm_inv
,   vector<float>& lpd
,   vector<float>& hpd
,   vector<float>& lpr
,   vector<float>& hpr
){
	//! Initialization for convenience
	vector<unsigned> row_ind(0);
    ind_initialize(row_ind, 0, h - prms.k, prms.p);

	vector<unsigned> column_ind(0);
    ind_initialize(column_ind, 0, w - prms.k, prms.p);

	const unsigned kHard_2 = prms.k * prms.k;
	vector<float> group_3D_table(d * kHard_2 * prms.N * column_ind.size());
	vector<float> wx_r_table;
	wx_r_table.reserve(d * column_ind.size());
	vector<float> hadamard_tmp(prms.N);

	//! For aggregation part
    memset(denominator, 0, sizeof(*denominator)*w*h*d);
    memset(basic_estimate, 0, sizeof(*basic_estimate)*w*h*d);

	//! Precompute Bloc-Matching
	vector<vector<pair<unsigned, unsigned> > > patch_table(column_ind.size(), std::vector<pair<unsigned,unsigned> >(prms.N));
	vector<unsigned> size_patch_table(column_ind.size());
	vector<float> distances(prms.N);

	vector<float> table_2D(prms.N * d * kHard_2, 0.0f);

    //! Loop on i_r
    for (unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++)
    {
        wx_r_table.clear();
        group_3D_table.clear();

        const unsigned i_r = row_ind[ind_i];

        //! Loop on j_r
        for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
        {
            //! Initialization
            const unsigned j_r = column_ind[ind_j];
            unsigned nSx_r = computeSimilarPatches(distances, patch_table[ind_j], index, d*(j_r + i_r*w), buffer, w, h, d, size_buffer, prms);
            size_patch_table[ind_j] = nSx_r;

            //! Update of table_2D
            if (prms.T_2D == DCT)
                dct_2d_process(table_2D, buffer, w, d, patch_table[ind_j], plan_2d,
                        prms.k, coef_norm);
            else if (prms.T_2D == BIOR)
                bior_2d_process(table_2D, buffer, w, d, patch_table[ind_j], 
                        prms.k, lpd, hpd);

            //! Build of the 3D group
            vector<float> group_3D(d * nSx_r * kHard_2, 0.0f);
            for(unsigned c = 0; c < d; c++)
                for (unsigned n = 0; n < nSx_r; n++)
                    for (unsigned k = 0; k < kHard_2; k++)
                        group_3D[n + k * nSx_r + c * kHard_2 * nSx_r] =
                            table_2D[k + n * kHard_2 + c * kHard_2 * prms.N];

            //! HT filtering of the 3D group
            vector<float> weight_table(d);
            if(prms.T_3D == HADAMARD)
                ht_filtering_hadamard(group_3D, hadamard_tmp, nSx_r, prms.k, d, sigma, prms.lambda3D, weight_table);
            else
                ht_filtering_haar(group_3D, hadamard_tmp, nSx_r, prms.k, d, sigma, prms.lambda3D, weight_table);

            //! Save the 3D group. The DCT 2D inverse will be done after.
            for (unsigned c = 0; c < d; c++)
                for (unsigned n = 0; n < nSx_r; n++)
                    for (unsigned k = 0; k < kHard_2; k++)
                        group_3D_table.push_back(group_3D[n + k * nSx_r + c * kHard_2 * nSx_r]);

            //! Save weighting
            for (unsigned c = 0; c < d; c++)
                wx_r_table.push_back(weight_table[c]);

        } //! End of loop on j_r

        //!  Apply 2D inverse transform
        if (prms.T_2D == DCT)
            dct_2d_inv(group_3D_table, prms.k, prms.N * d * column_ind.size(),
                    coef_norm_inv, plan_2d_inv);
        else if (prms.T_2D == BIOR)
            bior_2d_inv(group_3D_table, prms.k, lpr, hpr);


        //! Registration of the weighted estimation
        unsigned dec = 0;
        for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
        {
            const unsigned nSx_r = size_patch_table[ind_j];
            for (unsigned c = 0; c < d; c++)
            {
                for (unsigned n = 0; n < nSx_r; n++)
                {
                    //! Only aggregates patches from the current frame
                    if(patch_table[ind_j][n].first == index)
                    {
                        for (unsigned p = 0; p < prms.k; p++)
                            for (unsigned q = 0; q < prms.k; q++)
                            {
                                basic_estimate[patch_table[ind_j][n].second + d*(p + q*w)] += kaiser_window[p * prms.k + q] * wx_r_table[c + ind_j * d]
                                    * group_3D_table[p * prms.k + q + n * kHard_2 + c * kHard_2 * nSx_r + dec];
                                denominator[patch_table[ind_j][n].second + d*(p + q*w)] += kaiser_window[p * prms.k + q] * wx_r_table[c + ind_j * d];
                            }
                    }
                }
            }

            dec += nSx_r * d * kHard_2;
        }

    } //! End of loop on i_r
}

/**
 * @brief Run the final process of BM3D (2nd step). The result
 *        is contained in final_estimate. 
 *        
 * @param sigma: value of assumed noise of the video to denoise;
 * @param buffer_input: vector containing the noisy frames;
 * @param buffer_basic: vector containing the basic frames;
 * @param final_estimate: will be the denoised frame;
 * @param denominator: will containes the aggregation coefficients;
 * @param w: width of a frame;
 * @param h: height of a frame;
 * @param d: number of channels of a frame;
 * @param prms: parameters;
 * @param plan_2d: fftw plan for the dct transform;
 * @param plan_2d_inv: fftw plan for the inverse dct transform;
 * @param index: index of the frame in the buffer;
 * @param size_buffer: current size of the buffer;
 * @param kaiser_window: coefficient for the kaiser window;
 * @param coef_norm: normalization coeffs;
 * @param coef_norm_inv: inverse normalization coeffs;
 * @param lpd, hpd, lpr, hpr: bior basis;
 *
 * @return none.
 **/
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
,   vector<float>& kaiser_window
,   vector<float>& coef_norm
,   vector<float>& coef_norm_inv
,   vector<float>& lpd
,   vector<float>& hpd
,   vector<float>& lpr
,   vector<float>& hpr
){
	//! Initialization for convenience
	vector<unsigned> row_ind(0);
		ind_initialize(row_ind, 0, (h - prms.k), prms.p);

	vector<unsigned> column_ind(0);
		ind_initialize(column_ind, 0, (w - prms.k), prms.p);

	const unsigned kWien_2 = prms.k * prms.k;
	vector<float> group_3D_table(d * kWien_2 * prms.N * column_ind.size());
	vector<float> wx_r_table;
	wx_r_table.reserve(d * column_ind.size());
	vector<float> tmp(prms.N);

	//! For aggregation part
    memset(denominator, 0, sizeof(*denominator)*w*h*d);
    memset(final_estimate, 0, sizeof(*final_estimate)*w*h*d);

	//! Precompute Bloc-Matching
	vector<vector<pair<unsigned,unsigned> > > patch_table(column_ind.size(), std::vector<pair<unsigned,unsigned> >(prms.N));
	vector<unsigned> size_patch_table(column_ind.size());
	vector<float> distances(prms.N);

	//! DCT_table_2D[p * N + q + (i * width + j) * kWien_2 + c * (2 * ns + 1) * width * kWien_2]
	vector<float> table_2D_vid(prms.N * d * kWien_2, 0.0f);
	vector<float> table_2D_est(prms.N * d * kWien_2, 0.0f);

    //! Loop on i_r
    for (unsigned ind_i = 0; ind_i < row_ind.size(); ind_i++)
    {
        wx_r_table.clear();
        group_3D_table.clear();

        const unsigned i_r = row_ind[ind_i];

        //! Loop on j_r
        for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
        {
            const unsigned j_r = column_ind[ind_j];
            //! Number of similar patches
            unsigned nSx_r = computeSimilarPatches(distances, patch_table[ind_j], index, d*(j_r + i_r*w), buffer_basic, w, h, d, size_buffer, prms);
            size_patch_table[ind_j] = nSx_r;

            //! Update of DCT_table_2D
            if (prms.T_2D == DCT)
            {
                dct_2d_process(table_2D_vid, buffer_input, w, d, patch_table[ind_j], plan_2d,
                        prms.k, coef_norm);
                dct_2d_process(table_2D_est, buffer_basic, w, d, patch_table[ind_j], plan_2d,
                        prms.k, coef_norm);
            }
            else if (prms.T_2D == BIOR)
            {
                bior_2d_process(table_2D_vid, buffer_input, w, d, patch_table[ind_j], 
                        prms.k, lpd, hpd);
                bior_2d_process(table_2D_est, buffer_basic, w, d, patch_table[ind_j],
                        prms.k, lpd, hpd);
            }

            //! Build of the 3D group
            vector<float> group_3D_est(d * nSx_r * kWien_2, 0.0f);
            vector<float> group_3D_vid(d * nSx_r * kWien_2, 0.0f);
            for (unsigned c = 0; c < d; c++)
                for (unsigned n = 0; n < nSx_r; n++)
                {
                    for (unsigned k = 0; k < kWien_2; k++)
                    {
                        group_3D_est[n + k * nSx_r + c * kWien_2 * nSx_r] =
                            table_2D_est[k + n * kWien_2 + c * kWien_2 * prms.N];
                        group_3D_vid[n + k * nSx_r + c * kWien_2 * nSx_r] =
                            table_2D_vid[k + n * kWien_2 + c * kWien_2 * prms.N];
                    }
                }

            //! Wiener filtering of the 3D group
            vector<float> weight_table(d);
            if(prms.T_3D == HADAMARD)
                wiener_filtering_hadamard(group_3D_vid, group_3D_est, tmp, nSx_r, prms.k,
                        d, sigma, weight_table);
            else
                wiener_filtering_haar(group_3D_vid, group_3D_est, tmp, nSx_r, prms.k,
                        d, sigma, weight_table);

            //! Save the 3D group. The DCT 2D inverse will be done after.
            for (unsigned c = 0; c < d; c++)
                for (unsigned n = 0; n < nSx_r; n++)
                    for (unsigned k = 0; k < kWien_2; k++)
                        group_3D_table.push_back(group_3D_est[n + k * nSx_r + c * kWien_2 * nSx_r]);

            //! Save weighting
            for (unsigned c = 0; c < d; c++)
                wx_r_table.push_back(weight_table[c]);

        } //! End of loop on j_r

        //!  Apply 2D dct inverse
        if (prms.T_2D == DCT)
            dct_2d_inv(group_3D_table, prms.k, prms.N * d * column_ind.size(),
                    coef_norm_inv, plan_2d_inv);
        else if (prms.T_2D == BIOR)
            bior_2d_inv(group_3D_table, prms.k, lpr, hpr);

        //! Registration of the weighted estimation
        unsigned dec = 0;
        for (unsigned ind_j = 0; ind_j < column_ind.size(); ind_j++)
        {
            const unsigned nSx_r = size_patch_table[ind_j];
            for (unsigned c = 0; c < d; c++)
            {
                for (unsigned n = 0; n < nSx_r; n++)
                {
                    if(patch_table[ind_j][n].first == index)
                    {
                        for (unsigned p = 0; p < prms.k; p++)
                            for (unsigned q = 0; q < prms.k; q++)
                            {
                                final_estimate[patch_table[ind_j][n].second + d*(p*w + q)] += kaiser_window[p * prms.k + q] * wx_r_table[c + ind_j * d]
                                    * group_3D_table[p * prms.k + q + n * kWien_2 + c * kWien_2 * nSx_r + dec];
                                denominator[patch_table[ind_j][n].second + d*(p*w + q)] += kaiser_window[p * prms.k + q] * wx_r_table[c + ind_j * d];
                            }
                    }
                }
            }
            dec += nSx_r * d * kWien_2;
        }
	}
}

/**
 * @brief Computes the distance between to patches
 *
 * @param patch1, patch2 : indexes of the two patch to be compared;
 * @param frames: vector of frames;
 * @param w, d: width and number of channels of the frames;
 * @param patch_table: contains the indexes of the chosen patches;
 * @param size_buffer : size of the buffer of frames;
 * @param sizePatch : size of patches
 *
 * @return the distance of the patch
 **/
inline float patchDistance(
	unsigned patch1
, 	unsigned patch2
, 	vector<float*>& frames
,   const int w
,   const int d
,   const int size_buffer
, 	int sizePatch
){
	const int sPx = sizePatch;
    unsigned t1,p1;
    unsigned t2,p2;

    t1 = patch1 % size_buffer;
    t2 = patch2 % size_buffer;

    p1 = patch1 / size_buffer;
    p2 = patch2 / size_buffer;

	float dist = 0.f, dif;
	for (unsigned hc = 0; hc < d; ++hc)
    for (unsigned hy = 0; hy < sPx; hy++)
    for (unsigned hx = 0; hx < sPx; hx++)
        dist += (dif = (frames[t1][p1 + d*(hx + hy*w)] - frames[t2][p2 + d*(hx + hy*w)])) * dif;
	return dist / (sPx * sPx * d) / (255.f*255.f);
}

/**
 * @brief Computes the nearest patches in a given frame
 *
 * @param pidx : patch on which the search is centered;
 * @param rpidx : reference patch;
 * @param s: size of the search region;
 * @param k: size of the pathes;
 * @param Nb: number of patches to keep at the end;
 * @param d: bias for patches at the same spatial position than the reference;
 * @param frames: vector of frames;
 * @param w, h, d: width, height and number of channels of the frames;
 * @param size_buffer : size of the buffer of frames;
 * @param alreadySeen: index of the patch that has already been considered;
 * @param bestPatches: will contain the index of the best patches 
 *                      and their distance to the reference;
 **/
inline void localSearch(
	unsigned pidx
, 	unsigned rpidx
, 	unsigned s
, 	int k
, 	unsigned Nb
,	float d
, 	vector<float*>& frames
,   const int w
,   const int h
,   const int c
,   const int size_buffer
, 	std::unordered_map<unsigned, int>& alreadySeen
, 	std::vector<std::pair<float, unsigned> >& bestPatches
){
	int sWx = s;
	int sWy = s;
	const int sPx = k;

    unsigned rp;
    rp = rpidx / size_buffer;

    unsigned ct, cp;
    ct = pidx % size_buffer;
    cp = pidx / size_buffer;
    unsigned px, py;
    px = cp % w;
    py = cp / w;

	unsigned rangex[2];
	unsigned rangey[2];

	rangex[0] = std::max(0, (int)px - (sWx-1)/2);
	rangey[0] = std::max(0, (int)py - (sWy-1)/2);

	rangex[1] = std::min(w - sPx, (int)px + (sWx-1)/2);
	rangey[1] = std::min(h - sPx, (int)py + (sWy-1)/2);

	//! Redefine size of search range
	sWx = rangex[1] - rangex[0] + 1;
	sWy = rangey[1] - rangey[0] + 1;

	std::vector<std::pair<float, unsigned> > distance;
	distance.reserve(sWx*sWy);

	//! Compute distance between patches in search range
	for (unsigned qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
		for (unsigned qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
		{
            unsigned currentPatch = c*(qx + qy*w)*size_buffer + ct;

			//! Save distance and corresponding patch index
			int seen = (alreadySeen[currentPatch]++);
			if(seen == 0)
				distance.push_back((rp == (qx + qy*w)) ? std::make_pair(patchDistance(rpidx, currentPatch, frames, w, c, size_buffer, sPx) - d, currentPatch):std::make_pair(patchDistance(rpidx, currentPatch, frames, w, c, size_buffer, sPx), currentPatch));
		}

	int nbCandidates = std::min(Nb, (unsigned)distance.size());
	std::partial_sort(distance.begin(), distance.begin() + nbCandidates,
			distance.end(), comparaisonFirst);
	for(unsigned ix = 0; ix < nbCandidates; ++ix)
		bestPatches.push_back(distance[ix]);
}

/**
 * @brief Computes the nearest patches using the block matching principle of VBM3D
 *
 * @param output : will contain the distance to the best patches;
 * @param index : will contain the indexes to the best patches;
 * @param idx_curr_frame : index of the frame in which the search start;
 * @param pidx : spatial position of the reference patch
 * @param frames : vector of frames;
 * @param w, h, c : width, height and number of channels of the frames;
 * @param size_buffer : size of the buffer of frames;
 * @param prms : parameters
 **/
int computeSimilarPatches(
	std::vector<float>& output
,	std::vector<pair<unsigned,unsigned> >& index
,   int idx_curr_frame
,	unsigned pidx
, 	vector<float*>& frames
,   const int w
,   const int h
,   const int c
,   const int size_buffer
,	const Parameters& prms
){
	std::vector<unsigned> tempMatchesPre(prms.Nb);
	std::vector<std::pair<float, unsigned> > bestPatches;
	bestPatches.reserve(prms.Nb*(prms.Nf+1));
	std::vector<std::pair<float, unsigned> > frameBestPatches;
	frameBestPatches.reserve(prms.Nb*prms.Nb);
	std::unordered_map<unsigned, int> alreadySeen;

    int ridx = pidx * size_buffer + idx_curr_frame;

	//! Search in the current frame
	localSearch(ridx, ridx, prms.Ns, prms.k, prms.Nb, prms.d, frames, w, h, c, size_buffer, alreadySeen, bestPatches);
	for(unsigned ix = 0; ix < prms.Nb; ++ix)
		tempMatchesPre[ix] = bestPatches[ix].second;

	//! Search in the previous frames (centered on the matches)
	int finalFrame = (idx_curr_frame - min(size_buffer, (int) prms.Nf) + 1 + size_buffer) % size_buffer;
    int nextFrame = (idx_curr_frame + size_buffer - 1) % size_buffer;
    bool over = (finalFrame == idx_curr_frame) ? true : false;
    while(!over)
	{
		frameBestPatches.clear();
		for(unsigned currentTempMatch = 0; currentTempMatch < prms.Nb; ++currentTempMatch)
		{
            unsigned currentMatch = (tempMatchesPre[currentTempMatch] / size_buffer) * size_buffer + nextFrame;
            localSearch(currentMatch, ridx, prms.Npr, prms.k, prms.Nb, prms.d, frames, w, h, c, size_buffer, alreadySeen, frameBestPatches);
		}

		int nbCandidates = std::min(prms.Nb, (unsigned)frameBestPatches.size());
		std::partial_sort(frameBestPatches.begin(), frameBestPatches.begin() + nbCandidates,
				frameBestPatches.end(), comparaisonFirst);
		for(unsigned ix = 0; ix < nbCandidates; ++ix)
		{
			tempMatchesPre[ix] = frameBestPatches[ix].second;
			bestPatches.push_back(frameBestPatches[ix]);
		}
        if(nextFrame == finalFrame)
            over = true;
        else
            nextFrame = (nextFrame + size_buffer - 1) % size_buffer;
	}

	const unsigned nSimP = std::min(prms.N, (unsigned)bestPatches.size());

	std::partial_sort(bestPatches.begin(), bestPatches.begin() + nSimP,
			bestPatches.end(), comparaisonFirst);

	for (unsigned n = 0; n < nSimP; n++)
	{
		output[n] = bestPatches[n].first;
		index[n] = make_pair(bestPatches[n].second % size_buffer, bestPatches[n].second / size_buffer);
	}

	unsigned ind_thresh = nSimP - 1;
	while((output[ind_thresh] > prms.tau) && (ind_thresh > 0))
		ind_thresh--;

    int candidates = closest_power_of_2(ind_thresh+1);

#ifdef MTRICK
    // Artificially adds a candidate when there's only the reference patch left 
    if(candidates == 1)
    {
        candidates = 2;
        output[1] = output[0];
        index[1] = index[0];
    }
#endif

	return candidates;
}

/**
 * @brief Precompute a 2D DCT transform on all patches contained in
 *        a part of the video.
 *
 * @param DCT_table_2D : will contain the 2d DCT transform for all
 *        chosen patches;
 * @param buffer : vector containing the frames;
 * @param w, d: width and number of channels of the frames;
 * @param patch_table: contains the indexes of the chosen patches;
 * @param plan : for convenience. Used by fftw;
 * @param kHW : size of patches (kHW x kHW). MUST BE A POWER OF 2 !!!
 * @param lpd : low pass filter of the forward bior1.5 2d transform;
 * @param hpd : high pass filter of the forward bior1.5 2d transform.
 **/
void dct_2d_process(
    vector<float> &DCT_table_2D
,   vector<float*> const& buffer
,   const int w
,   const int d
,   vector<pair<unsigned,unsigned> > const& patch_table 
,   fftwf_plan * plan
,   const unsigned kHW
,   vector<float> const& coef_norm
){
	//! Declarations
	const unsigned kHW_2 = kHW * kHW;
	const unsigned size = d * kHW_2 * patch_table.size();

	//! Allocating Memory
	float* vec = (float*) fftwf_malloc(size * sizeof(float));
	float* dct = (float*) fftwf_malloc(size * sizeof(float));

	for (unsigned c = 0; c < d; c++)
	{
		const unsigned dc_p = c * kHW_2 * patch_table.size();
		for(unsigned n = 0; n < patch_table.size(); ++n)
		{
            for (unsigned p = 0; p < kHW; p++)
                for (unsigned q = 0; q < kHW; q++)
                    vec[p * kHW + q + dc_p + n * kHW_2] =
                        buffer[patch_table[n].first][patch_table[n].second + d*(p*w+q)];
		}
	}

	//! Process of all DCTs
	fftwf_execute_r2r(*plan, vec, dct);
	fftwf_free(vec);

	//! Getting the result
	for (unsigned c = 0; c < d; c++)
	{
		const unsigned dc_p = c * kHW_2 * patch_table.size();
		for(unsigned n = 0; n < patch_table.size(); ++n)
            for (unsigned k = 0; k < kHW_2; k++)
                DCT_table_2D[dc_p + n * kHW_2 + k] =
                    dct[dc_p + n * kHW_2 + k] * coef_norm[k];
	}
	fftwf_free(dct);
}

/**
 * @brief Precompute a 2D bior1.5 transform on all patches contained in
 *        a part of the video.
 *
 * @param bior_table_2D : will contain the 2d bior1.5 transform for all
 *        chosen patches;
 * @param buffer : vector containing the frames;
 * @param w, d: width and number of channels of the frames;
 * @param patch_table: contains the indexes of the chosen patches;
 * @param kHW : size of patches (kHW x kHW). MUST BE A POWER OF 2 !!!
 * @param lpd : low pass filter of the forward bior1.5 2d transform;
 * @param hpd : high pass filter of the forward bior1.5 2d transform.
 **/
void bior_2d_process(
    vector<float> &bior_table_2D
,   vector<float*> const& buffer
,   const int w
,   const int d
,   vector<pair<unsigned,unsigned> > const& patch_table
,   const unsigned kHW
,   vector<float> &lpd
,   vector<float> &hpd
){
	//! Declarations
	const unsigned kHW_2 = kHW * kHW;

	//! If i_r == ns, then we have to process all Bior1.5 transforms
	for (unsigned c = 0; c < d; c++)
	{
		const unsigned dc_p = c * kHW_2 * patch_table.size();
		for(unsigned n = 0; n < patch_table.size(); ++n)
            bior_2d_forward(buffer, w, bior_table_2D, kHW, patch_table[n].second + c, patch_table[n].first, d, dc_p + n * kHW_2, lpd, hpd);
	}
}

/**
 * @brief HT filtering using Welsh-Hadamard transform (do only third
 *        dimension transform, Hard Thresholding and inverse transform).
 *
 * @param group_3D : contains the 3D block for a reference patch;
 * @param tmp: allocated vector used in Hadamard transform for convenience;
 * @param nSx_r : number of similar patches to a reference one;
 * @param kHW : size of patches (kHW x kHW);
 * @param chnls : number of channels of the video;
 * @param sigma : value of noise;
 * @param lambdaHard3D : value of thresholding;
 * @param weight_table: the weighting of this 3D group for each channel;
 *
 * @return none.
 **/
void ht_filtering_hadamard(
    vector<float> &group_3D
,   vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kHard
,   const unsigned chnls
,   float sigma
,   const float lambdaHard3D
,   vector<float> &weight_table
){
	//! Declarations
	const unsigned kHard_2 = kHard * kHard;
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 0.0f;
	const float coef_norm = sqrtf((float) nSx_r);
	const float coef = 1.0f / (float) nSx_r;

	//! Process the Welsh-Hadamard transform on the 3rd dimension
	for (unsigned n = 0; n < kHard_2 * chnls; n++)
		hadamard_transform(group_3D, tmp, nSx_r, n * nSx_r);

	//! Hard Thresholding
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc = c * nSx_r * kHard_2;
		const float T = lambdaHard3D * sigma * coef_norm;
		for (unsigned k = 0; k < kHard_2 * nSx_r; k++)
		{
#ifdef DCTHRESH
            if (fabs(group_3D[k + dc]) > T)
#else
            if (k < nSx_r || fabs(group_3D[k + dc]) > T)
#endif
				weight_table[c]++;
			else
				group_3D[k + dc] = 0.0f;
		}
	}

	//! Process of the Welsh-Hadamard inverse transform
	for (unsigned n = 0; n < kHard_2 * chnls; n++)
		hadamard_transform(group_3D, tmp, nSx_r, n * nSx_r);

	for (unsigned k = 0; k < group_3D.size(); k++)
		group_3D[k] *= coef;

	//! Weight for aggregation
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 1.0f / (float) (sigma * sigma * weight_table[c]);
}

/**
 * @brief HT filtering using Haar (do only third
 *        dimension transform, Hard Thresholding and inverse transform).
 *
 * @param group_3D : contains the 3D block for a reference patch;
 * @param tmp: allocated vector used in Hadamard transform for convenience;
 * @param nSx_r : number of similar patches to a reference one;
 * @param kHW : size of patches (kHW x kHW);
 * @param chnls : number of channels of the video;
 * @param sigma : value of noise;
 * @param lambdaHard3D : value of thresholding;
 * @param weight_table: the weighting of this 3D group for each channel;
 *        otherwise.
 *
 * @return none.
 **/
void ht_filtering_haar(
    vector<float> &group_3D
,   vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kHard
,   const unsigned chnls
,   float sigma
,   const float lambdaHard3D
,   vector<float> &weight_table
){
	//! Declarations
	const unsigned kHard_2 = kHard * kHard;
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 0.0f;

	//! Process the Haar transform on the 3rd dimension
	for (unsigned n = 0; n < kHard_2 * chnls; n++)
		haar_forward(group_3D, tmp, nSx_r, n * nSx_r);

	//! Hard Thresholding
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc = c * nSx_r * kHard_2;
		const float T = lambdaHard3D * sigma;
		for (unsigned k = 0; k < kHard_2 * nSx_r; k++)
		{
#ifdef DCTHRESH
            if (fabs(group_3D[k + dc]) > T)
#else
            if (k < nSx_r || fabs(group_3D[k + dc]) > T)
#endif
				weight_table[c]++;
			else
				group_3D[k + dc] = 0.0f;
		}
	}

	//! Process of the Haar inverse transform
	for (unsigned n = 0; n < kHard_2 * chnls; n++)
		haar_inverse(group_3D, tmp, nSx_r, n * nSx_r);

	//! Weight for aggregation
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 1.f / (float) (sigma * sigma * weight_table[c]);
}

/**
 * @brief Wiener filtering using Hadamard transform.
 *
 * @param group_3D_vid : contains the 3D block built on the noisy sequence;
 * @param group_3D_est : contains the 3D block built on the basic sequence;
 * @param tmp: allocated vector used in hadamard transform for convenience;
 * @param nSx_r : number of similar patches to a reference one;
 * @param kWien : size of patches (kWien x kWien);
 * @param chnls : number of channels of the video;
 * @param sigma : value of noise;
 * @param weight_table: the weighting of this 3D group for each channel;
 *
 * @return none.
 **/
void wiener_filtering_hadamard(
    vector<float> &group_3D_vid
,   vector<float> &group_3D_est
,   vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kWien
,   const unsigned chnls
,   float sigma
,   vector<float> &weight_table
){
	//! Declarations
	const unsigned kWien_2 = kWien * kWien;
	const float coef = 1.0f / (float) nSx_r;

	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 0.0f;

	//! Process the Welsh-Hadamard transform on the 3rd dimension
	for (unsigned n = 0; n < kWien_2 * chnls; n++)
	{
		hadamard_transform(group_3D_vid, tmp, nSx_r, n * nSx_r);
		hadamard_transform(group_3D_est, tmp, nSx_r, n * nSx_r);
	}

	//! Wiener Filtering
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc = c * nSx_r * kWien_2;
#ifdef DCWIENER
		for (unsigned k = 0; k < kWien_2 * nSx_r; k++)
#else
        for (unsigned k = nSx_r; k < kWien_2 * nSx_r; k++)
#endif
		{
			float value = group_3D_est[dc + k] * group_3D_est[dc + k] * coef;
			value /= (value + sigma * sigma);
			group_3D_est[k + dc] = group_3D_vid[k + dc] * value * coef;
			weight_table[c] += (value*value);
		}
#ifndef DCWIENER
        // Add the weight corresponding to the DC components that was not thresholded
        weight_table[c] += nSx_r; 
#endif
	}

	//! Process of the Welsh-Hadamard inverse transform
	for (unsigned n = 0; n < kWien_2 * chnls; n++)
		hadamard_transform(group_3D_est, tmp, nSx_r, n * nSx_r);

	//! Weight for aggregation
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = (weight_table[c] > 0.0f ? 1.0f / (float)
				(sigma * sigma * weight_table[c]) : 1.0f);
}


/**
 * @brief Wiener filtering using Haar transform.
 *
 * @param group_3D_vid : contains the 3D block built on the noisy sequence;
 * @param group_3D_est : contains the 3D block built on the basic sequence;
 * @param tmp: allocated vector used in hadamard transform for convenience;
 * @param nSx_r : number of similar patches to a reference one;
 * @param kWien : size of patches (kWien x kWien);
 * @param chnls : number of channels of the video;
 * @param sigma : value of noise;
 * @param weight_table: the weighting of this 3D group for each channel;
 *
 * @return none.
 **/
void wiener_filtering_haar(
    vector<float> &group_3D_vid
,   vector<float> &group_3D_est
,   vector<float> &tmp
,   const unsigned nSx_r
,   const unsigned kWien
,   const unsigned chnls
,   float sigma
,   vector<float> &weight_table
){
	//! Declarations
	const unsigned kWien_2 = kWien * kWien;

	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 0.0f;

	//! Process the Haar transform on the 3rd dimension
	for (unsigned n = 0; n < kWien_2 * chnls; n++)
	{
		haar_forward(group_3D_vid, tmp, nSx_r, n * nSx_r);
		haar_forward(group_3D_est, tmp, nSx_r, n * nSx_r);
	}

	//! Wiener Filtering
	for (unsigned c = 0; c < chnls; c++)
	{
		const unsigned dc = c * nSx_r * kWien_2;
#ifdef DCWIENER
		for (unsigned k = 0; k < kWien_2 * nSx_r; k++)
#else
        for (unsigned k = nSx_r; k < kWien_2 * nSx_r; k++)
#endif
		{
			float value = group_3D_est[dc + k] * group_3D_est[dc + k];
			value /= (value + sigma * sigma);
			group_3D_est[k + dc] = group_3D_vid[k + dc] * value;
			weight_table[c] += (value*value);
		}
#ifndef DCWIENER
        // Add the weight corresponding to the DC components that were not thresholded
        weight_table[c] += nSx_r; 
#endif
	}

	//! Process of the Welsh-Hadamard inverse transform
	for (unsigned n = 0; n < kWien_2 * chnls; n++)
		haar_inverse(group_3D_est, tmp, nSx_r, n * nSx_r);

	//! Weight for aggregation
	for (unsigned c = 0; c < chnls; c++)
		weight_table[c] = 1.f / (float) (sigma * sigma * weight_table[c]);
}

/**
 * @brief Apply 2D dct inverse to a lot of patches.
 *
 * @param group_3D_table: contains a huge number of patches;
 * @param kHW : size of patch;
 * @param coef_norm_inv: contains normalization coefficients;
 * @param plan : for convenience. Used by fftw.
 *
 * @return none.
 **/
void dct_2d_inv(
		vector<float> &group_3D_table
		,   const unsigned kHW
		,   const unsigned N
		,   vector<float> const& coef_norm_inv
		,   fftwf_plan * plan
		){
	//! Declarations
	const unsigned kHW_2 = kHW * kHW;
	const unsigned size = kHW_2 * N;
	const unsigned Ns   = group_3D_table.size() / kHW_2;

	//! Allocate Memory
	float* vec = (float*) fftwf_malloc(size * sizeof(float));
	float* dct = (float*) fftwf_malloc(size * sizeof(float));

	//! Normalization
	for (unsigned n = 0; n < Ns; n++)
        for (unsigned k = 0; k < kHW_2; k++)
            dct[k + n * kHW_2] = group_3D_table[k + n * kHW_2] * coef_norm_inv[k];

	//! 2D dct inverse
	fftwf_execute_r2r(*plan, dct, vec);
	fftwf_free(dct);

	//! Getting the result + normalization
	const float coef = 1.0f / (float)(kHW * 2);
	for (unsigned k = 0; k < group_3D_table.size(); k++)
		group_3D_table[k] = coef * vec[k];

	//! Free Memory
	fftwf_free(vec);
}

void bior_2d_inv(
		vector<float> &group_3D_table
		,   const unsigned kHW
		,   vector<float> const& lpr
		,   vector<float> const& hpr
		){
	//! Declarations
	const unsigned kHW_2 = kHW * kHW;
	const unsigned N = group_3D_table.size() / kHW_2;

	//! Bior process
	for (unsigned n = 0; n < N; n++)
        bior_2d_inverse(group_3D_table, kHW, n * kHW_2, lpr, hpr);
}

/** ----------------- **/
/** - Preprocessing - **/
/** ----------------- **/
/**
 * @brief Preprocess
 *
 * @param kaiser_window[kHW * kHW]: Will contain values of a Kaiser Window;
 * @param coef_norm: Will contain values used to normalize the 2D DCT;
 * @param coef_norm_inv: Will contain values used to normalize the 2D DCT;
 * @param kHW: size of patches (need to be 8 or 12).
 *
 * @return none.
 **/
void preProcess(
		vector<float> &kaiserWindow
		,   vector<float> &coef_norm
		,   vector<float> &coef_norm_inv
		,   const unsigned kHW
	       ){
	//! Kaiser Window coefficients
	if(kHW == 4)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.4055f;
		kaiserWindow[1 + kHW * 0] = 0.4055f; kaiserWindow[1 + kHW * 1] = 0.8544f;

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i < kHW / 2; i++)
			for (unsigned j = kHW / 2; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else if (kHW == 6)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.3368f; kaiserWindow[0 + kHW * 2] = 0.4265f;
		kaiserWindow[1 + kHW * 0] = 0.3368f; kaiserWindow[1 + kHW * 1] = 0.5893f; kaiserWindow[1 + kHW * 2] = 0.7464f;
		kaiserWindow[2 + kHW * 0] = 0.4265f; kaiserWindow[2 + kHW * 1] = 0.7464f; kaiserWindow[2 + kHW * 2] = 0.9454f;

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i < kHW / 2; i++)
			for (unsigned j = kHW / 2; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else if (kHW == 7)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.3151f; kaiserWindow[0 + kHW * 2] = 0.4055f; kaiserWindow[0 + kHW * 3] = 0.4387f;
		kaiserWindow[1 + kHW * 0] = 0.3151f; kaiserWindow[1 + kHW * 1] = 0.5161f; kaiserWindow[1 + kHW * 2] = 0.6640f; kaiserWindow[1 + kHW * 3] = 0.7184f;
		kaiserWindow[2 + kHW * 0] = 0.4055f; kaiserWindow[2 + kHW * 1] = 0.6640f; kaiserWindow[2 + kHW * 2] = 0.8544f; kaiserWindow[2 + kHW * 3] = 0.9243f;
		kaiserWindow[3 + kHW * 0] = 0.4387f; kaiserWindow[3 + kHW * 1] = 0.7184f; kaiserWindow[3 + kHW * 2] = 0.9243f; kaiserWindow[3 + kHW * 3] = 1.0000f; 

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i <= kHW / 2; i++)
			for (unsigned j = kHW / 2 + 1; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2 + 1; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else if (kHW == 8)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.2989f; kaiserWindow[0 + kHW * 2] = 0.3846f; kaiserWindow[0 + kHW * 3] = 0.4325f;
		kaiserWindow[1 + kHW * 0] = 0.2989f; kaiserWindow[1 + kHW * 1] = 0.4642f; kaiserWindow[1 + kHW * 2] = 0.5974f; kaiserWindow[1 + kHW * 3] = 0.6717f;
		kaiserWindow[2 + kHW * 0] = 0.3846f; kaiserWindow[2 + kHW * 1] = 0.5974f; kaiserWindow[2 + kHW * 2] = 0.7688f; kaiserWindow[2 + kHW * 3] = 0.8644f;
		kaiserWindow[3 + kHW * 0] = 0.4325f; kaiserWindow[3 + kHW * 1] = 0.6717f; kaiserWindow[3 + kHW * 2] = 0.8644f; kaiserWindow[3 + kHW * 3] = 0.9718f;

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i < kHW / 2; i++)
			for (unsigned j = kHW / 2; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else if (kHW == 12)
	{
		//! First quarter of the matrix
		kaiserWindow[0 + kHW * 0] = 0.1924f; kaiserWindow[0 + kHW * 1] = 0.2615f; kaiserWindow[0 + kHW * 2] = 0.3251f; kaiserWindow[0 + kHW * 3] = 0.3782f;  kaiserWindow[0 + kHW * 4] = 0.4163f;  kaiserWindow[0 + kHW * 5] = 0.4362f;
		kaiserWindow[1 + kHW * 0] = 0.2615f; kaiserWindow[1 + kHW * 1] = 0.3554f; kaiserWindow[1 + kHW * 2] = 0.4419f; kaiserWindow[1 + kHW * 3] = 0.5139f;  kaiserWindow[1 + kHW * 4] = 0.5657f;  kaiserWindow[1 + kHW * 5] = 0.5927f;
		kaiserWindow[2 + kHW * 0] = 0.3251f; kaiserWindow[2 + kHW * 1] = 0.4419f; kaiserWindow[2 + kHW * 2] = 0.5494f; kaiserWindow[2 + kHW * 3] = 0.6390f;  kaiserWindow[2 + kHW * 4] = 0.7033f;  kaiserWindow[2 + kHW * 5] = 0.7369f;
		kaiserWindow[3 + kHW * 0] = 0.3782f; kaiserWindow[3 + kHW * 1] = 0.5139f; kaiserWindow[3 + kHW * 2] = 0.6390f; kaiserWindow[3 + kHW * 3] = 0.7433f;  kaiserWindow[3 + kHW * 4] = 0.8181f;  kaiserWindow[3 + kHW * 5] = 0.8572f;
		kaiserWindow[4 + kHW * 0] = 0.4163f; kaiserWindow[4 + kHW * 1] = 0.5657f; kaiserWindow[4 + kHW * 2] = 0.7033f; kaiserWindow[4 + kHW * 3] = 0.8181f;  kaiserWindow[4 + kHW * 4] = 0.9005f;  kaiserWindow[4 + kHW * 5] = 0.9435f;
		kaiserWindow[5 + kHW * 0] = 0.4362f; kaiserWindow[5 + kHW * 1] = 0.5927f; kaiserWindow[5 + kHW * 2] = 0.7369f; kaiserWindow[5 + kHW * 3] = 0.8572f;  kaiserWindow[5 + kHW * 4] = 0.9435f;  kaiserWindow[5 + kHW * 5] = 0.9885f;

		//! Completing the rest of the matrix by symmetry
		for(unsigned i = 0; i < kHW / 2; i++)
			for (unsigned j = kHW / 2; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[i + kHW * (kHW - j - 1)];

		for (unsigned i = kHW / 2; i < kHW; i++)
			for (unsigned j = 0; j < kHW; j++)
				kaiserWindow[i + kHW * j] = kaiserWindow[kHW - i - 1 + kHW * j];
	}
	else
		for (unsigned k = 0; k < kHW * kHW; k++)
			kaiserWindow[k] = 1.0f;

	//! Coefficient of normalization for DCT II and DCT II inverse
	const float coef = 0.5f / ((float) (kHW));
	for (unsigned i = 0; i < kHW; i++)
		for (unsigned j = 0; j < kHW; j++)
		{
			if (i == 0 && j == 0)
			{
				coef_norm    [i * kHW + j] = 0.5f * coef;
				coef_norm_inv[i * kHW + j] = 2.0f;
			}
			else if (i * j == 0)
			{
				coef_norm    [i * kHW + j] = SQRT2_INV * coef;
				coef_norm_inv[i * kHW + j] = SQRT2;
			}
			else
			{
				coef_norm    [i * kHW + j] = 1.0f * coef;
				coef_norm_inv[i * kHW + j] = 1.0f;
			}
		}
}
