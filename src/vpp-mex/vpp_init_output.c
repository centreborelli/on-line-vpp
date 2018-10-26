#include <stdint.h>
#include "mex.h"

#include <stdio.h>
#include "vpp.h"

// FILE* vpp_init_output(const char* filename, int w, int h, int d);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) {
        mexErrMsgTxt("want only two arguments");
    }

    const char* filename = mxArrayToString(prhs[0]);

    int* sizes;
    {
        mxArray *lhs[1];
        mexCallMATLAB(1, lhs, 1, (void**)&(prhs[1]), "int32");
        sizes = (int*)mxGetData(lhs[0]);
    }

    int w = sizes[1];
    int h = sizes[0];
    int d = sizes[2];

    FILE* file = vpp_init_output(filename, w, h, d);
    mxFree((void*) filename);

    if (file) {
        mwSize dims[] = {4};
        plhs[0] = mxCreateNumericArray(1, dims, mxUINT64_CLASS, mxREAL);
        UINT64_T* pr = (UINT64_T *) mxGetData(plhs[0]);
        pr[1] = w;
        pr[0] = h;
        pr[2] = d;
        pr[3] = (UINT64_T) file;
    } else {
        plhs[0] = mxCreateLogicalMatrix(1, 1);
        *mxGetLogicals(plhs[0]) = 0;
    }
}


