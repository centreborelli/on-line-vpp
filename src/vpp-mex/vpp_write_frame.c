#include <stdint.h>
#include "mex.h"

#include <stdio.h>
#include "vpp.h"

// int vpp_write_frame(FILE* out, float* frame, int w, int h, int d);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) {
        mexErrMsgTxt("want only two arguments");
    }

    UINT64_T* dat;
    {
        mxArray *lhs[1];
        mexCallMATLAB(1, lhs, 1, (void**)&(prhs[0]), "uint64");
        dat = (UINT64_T*)mxGetData(lhs[0]);
    }

    int w = dat[1];
    int h = dat[0];
    int d = dat[2];
    FILE* file = (FILE*) dat[3];

    mwSize ndim = mxGetNumberOfDimensions(prhs[1]);
    const mwSize* dims = mxGetDimensions(prhs[1]);
    if (ndim >= 1 && dims[0] != h) {
        mexErrMsgTxt("invalid h");
    }
    if (ndim >= 2 && dims[1] != w) {
        mexErrMsgTxt("invalid w");
    }
    if (ndim >= 3 && dims[2] != d) {
        mexErrMsgTxt("invalid d");
    }

    float* array = malloc(sizeof(float)*w*h*d);

    float* pr;
    {
        mxArray *lhs[1];
        mexCallMATLAB(1, lhs, 1, (void**)&(prhs[1]), "single");
        pr = (float*)mxGetData(lhs[0]);
    }

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            for (int c = 0; c < d; c++) {
                array[(y * w + x) * d + c] = pr[(w * c + x) * h + y];
            }
        }
    }

    int ret = vpp_write_frame(file, array, w, h, d);

    if (!ret) {
        plhs[0] = mxCreateLogicalMatrix(1, 1);
        *mxGetLogicals(plhs[0]) = 0;
    } else {
        plhs[0] = mxCreateLogicalMatrix(1, 1);
        *mxGetLogicals(plhs[0]) = 1;
    }

    free(array);
}


