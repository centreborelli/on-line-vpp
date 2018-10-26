#include <stdint.h>
#include "mex.h"

#include <stdio.h>
#include "vpp.h"

// int vpp_read_frame(FILE* in, float* frame, int w, int h, int d);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1) {
        mexErrMsgTxt("want only one argument");
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

    float* array = malloc(sizeof(float)*w*h*d);
    int ret = vpp_read_frame(file, array, w, h, d);

    if (!ret) {
        plhs[0] = mxCreateLogicalMatrix(1, 1);
        *mxGetLogicals(plhs[0]) = 0;
    } else {
        mwSize dims[] = {h, w, d};
        plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
        float* pr = (float*) mxGetData(plhs[0]);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < d; c++) {
                    pr[(w * c + x) * h + y] = array[(y * w + x) * d + c];
                }
            }
        }
    }

    free(array);
}

