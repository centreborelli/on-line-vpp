#include "mex.h"

#include <stdio.h>
#include "vpp.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1) {
        mexErrMsgTxt("want only one argument");
    }

    UINT64_T* dat = mxGetData(prhs[0]);
    FILE* file = (FILE*) dat[3];
    fclose(file);
}



