#include "argparse/argparse.h"   // command line parser

#include <stdlib.h>
#include <math.h>       // nans (used as boundary value by bicubic interp)
#include <fftw3.h>      // computes dct
#include <omp.h>
#include <stdbool.h>

#include "vpp.h"
#include "nlkalman.c"

// main funcion [[[1

// 'usage' message in the command line
static const char *const usages[] = {
	"nlkalman-bwd [options] [[--] args]",
	"nlkalman-bwd [options]",
	NULL,
};

int main(int argc, const char *argv[])
{
	// parse command line [[[2

	// command line parameters and their defaults
	const char *nisy_path = NULL; // input streams from pipeline
	const char *flow_path = NULL;
	const char *occl_path = NULL;
	const char *deno_path = NULL; // output stream
	int fframe = 0, lframe = -1;
	float sigma = 0.f;
	bool verbose = false;
	struct nlkalman_params prms;
	prms.patch_sz     = -1; // -1 means automatic value
	prms.search_sz    = -1;
	prms.dista_th     = -1.;
	prms.beta_x       = -1.;
	prms.beta_t       = -1.;
	prms.dista_lambda = -1.;
	prms.pixelwise = false;

	// configure command line parser
	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Algorithm options"),
		OPT_STRING ('i', "nisy"  , &nisy_path, "noisy input path (printf format)"),
		OPT_STRING ('o', "flow"  , &flow_path, "backward flow path (printf format)"),
		OPT_STRING ('k', "occl"  , &occl_path, "flow occlusions mask (printf format)"),
		OPT_STRING ('d', "deno"  , &deno_path, "denoised output path (printf format)"),
		OPT_INTEGER('f', "first" , &fframe, "first frame"),
		OPT_INTEGER('l', "last"  , &lframe , "last frame"),
		OPT_FLOAT  ('s', "sigma" , &sigma, "noise standard dev"),
		OPT_INTEGER('p', "patch" , &prms.patch_sz, "patch size"),
		OPT_INTEGER('w', "search", &prms.search_sz, "search region radius"),
		OPT_FLOAT  ( 0 , "dth"   , &prms.dista_th, "patch distance threshold"),
		OPT_FLOAT  ( 0 , "beta_x", &prms.beta_x, "noise multiplier in spatial filtering"),
		OPT_FLOAT  ( 0 , "beta_t", &prms.beta_t, "noise multiplier in kalman filtering"),
		OPT_FLOAT  ( 0 , "lambda", &prms.dista_lambda, "noisy patch weight in patch distance"),
		OPT_BOOLEAN( 0 , "pixel" , &prms.pixelwise, "toggle pixel-wise denoising"),
		OPT_GROUP("Program options"),
		OPT_BOOLEAN('v', "verbose", &verbose, "verbose output"),
		OPT_END(),
	};

	// parse command line
	struct argparse argparse;
	argparse_init(&argparse, options, usages, 0);
	argparse_describe(&argparse, "\nA video denoiser based on non-local means.", "");
	argc = argparse_parse(&argparse, argc, argv);

	// default value for noise-dependent params
	nlkalman_default_params(&prms, sigma);

	// print parameters
	if (verbose)
	{
		printf("parameters:\n");
		printf("\tnoise  %f\n", sigma);
		printf("\t%s-wise mode\n", prms.pixelwise ? "pixel" : "patch");
		printf("\tpatch     %d\n", prms.patch_sz);
		printf("\tsearch    %d\n", prms.search_sz);
		printf("\tdth       %g\n", prms.dista_th);
		printf("\tlambda    %g\n", prms.dista_lambda);
		printf("\tbeta_x    %g\n", prms.beta_x);
		printf("\tbeta_t    %g\n", prms.beta_t);
		printf("\n");
#ifdef WEIGHTED_AGGREGATION
		printf("\tWEIGHTED_AGGREGATION ON\n");
#endif
	}

	// process first frame [[[2
	int w,h,c;
	FILE* nisy_in = vpp_init_input(nisy_path, &w, &h, &c);
	if (!nisy_in)
	{
		fprintf(stderr, "kalman: cannot initialize input '%s'\n", nisy_path);
		return EXIT_FAILURE;
	}

	FILE* deno_out = vpp_init_output(deno_path, w, h, c);
	if (!deno_out)
	{
		fprintf(stderr, "kalman: cannot initialize output '%s'\n", deno_path);
		return EXIT_FAILURE;
	}

	// allocate input and output frame buffers
	const int whc = w*h*c, wh2 = w*h*2;
	float* nisy1 = malloc(whc*sizeof(float));
	float* deno1 = malloc(whc*sizeof(float));

	// run denoising
	if (verbose) printf("processing frame %d\n", 0);
	vpp_read_frame(nisy_in, nisy1, w, h, c);
	nlkalman_frame(deno1, nisy1, NULL, w, h, c, sigma, prms, 0);
	vpp_write_frame(deno_out, deno1, w, h, c);

	// process following frames [[[2

	// load optical flow [[[3
	FILE* flow_in = NULL;
	if (flow_path)
	{
		if (verbose) printf("loading flow %s\n", flow_path);
		int w1, h1, c1;
		flow_in = vpp_init_input(flow_path, &w1, &h1, &c1);

		if (!flow_in)
		{
			fprintf(stderr, "kalman: cannot initialize input '%s'\n", flow_path);
			return free(nisy1), free(deno1), EXIT_FAILURE;
		}

		if (w*h != w1*h1 || c1 != 2)
		{
			fprintf(stderr, "Video and optical flow size missmatch\n");
			return free(nisy1), free(deno1), EXIT_FAILURE;
		}
	}

	// load occlusion masks [[[3
	FILE* occl_in = NULL;
	if (flow_path && occl_path)
	{
		if (verbose) printf("loading occl. mask %s\n", occl_path);
		int w1, h1, c1;
		occl_in = vpp_init_input(occl_path, &w1, &h1, &c1);

		if (!occl_in)
		{
			fprintf(stderr, "kalman: cannot initialize input '%s'\n", occl_path);
			return free(nisy1), free(deno1), EXIT_FAILURE;
		}

		if (w*h != w1*h1 || c1 != 1)
		{
			fprintf(stderr, "Video and occl. masks size missmatch\n");
			return free(nisy1), free(deno1), EXIT_FAILURE;
		}
	}


	// run denoiser [[[3
	float* warp0 = malloc(whc*sizeof(float));
	float* flow10 = (flow_in) ? malloc(wh2*sizeof(float)) : NULL;
	float* occl10 = (occl_in) ? malloc(w*h*sizeof(float)) : NULL;
	int f = 1;
	while (vpp_read_frame(nisy_in, nisy1, w, h, c))
	{
		if (verbose) printf("processing frame %d\n", f);

		// warp previous denoised frame
		if (flow10)
		{
			vpp_read_frame(flow_in, flow10, w, h, 2);
			if (occl10) vpp_read_frame(occl_in, occl10, w, h, 1);
			warp_bicubic_inplace(warp0, deno1, flow10, occl10, w, h, c);
		}
		else
			// copy without warping
			memcpy(warp0, deno1, whc*sizeof(float));

		// run denoising
		nlkalman_frame(deno1, nisy1, warp0, w, h, c, sigma, prms, f++);

		// write output to pipeline
		vpp_write_frame(deno_out, deno1, w, h, c);
	}

	// free mem and return [[[2
	free(deno1);
	free(nisy1);
	free(warp0);
	if (flow10) free(flow10);
	if (occl10) free(occl10);

	return EXIT_SUCCESS; // ]]]2
}

// vim:set foldmethod=marker:
// vim:set foldmarker=[[[,]]]:
