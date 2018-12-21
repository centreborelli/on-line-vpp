#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

/**
 * @brief Structures of parameters
 *
 **/

struct Parameters
{
	/// Type of the 2D tranform
	unsigned T_2D;
	/// Type of the 1D tranform 
	unsigned T_3D;
	/// Maximum number of similar patches
	unsigned Nmax;
	/// Number of frames forward (and backward) used during the search
	unsigned Rf;
	/// Size of the search region in the reference frame
	unsigned Rr;
	/// Size of the search region in the other frame
	unsigned Rp;
	/// Maximum number of matches kept for a local search region
	unsigned Nl;
	/// Size of the patch (spatial)
	unsigned k;
	/// Step
	unsigned st;
	/// Correcting parameter in the distance computation
	float d;
	/// Threshold if it's a hard thresholding step
	float lambda;
	/// Distance threshold
	float tau;

	/// Border of the tile when using the multithreading
	int n = 16;
};

#endif
