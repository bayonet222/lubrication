#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int M, i;
	double *h, *H;
	
	M = mxGetM(prhs[0]);
	
	plhs[0] = mxCreateDoubleMatrix(M*2-1, 1, mxREAL);
	h = mxGetPr(plhs[0]);
	
	H = mxGetPr(prhs[0]);

	for (i=1; i<M-1; i++){
		h[i*2-1] = H[i];
		h[i*2-2] = 0.5*(H[i]+H[i-1]);
	}
	
	return;
}