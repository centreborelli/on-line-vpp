#include "mex.h"

#include <stdio.h>
#include "vpp.h"

// FILE* vpp_init_input(const char* filename, int* w, int* h, int* d);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1) {
        mexErrMsgTxt("want only one argument");
    }

    const char* filename = mxArrayToString(prhs[0]);

    int w, h, d;
    FILE* file = vpp_init_input(filename, &w, &h, &d);
    mxFree((void*) filename);

    if (file) {
        long dims[] = {4};
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

