#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int M, i;
	double *h, *H;
	
	M = mxGetM(prhs[0]);
	
	plhs[0] = mxCreateDoubleMatrix((M+1)/2, 1, mxREAL);
	H = mxGetPr(plhs[0]);
	
	h = mxGetPr(prhs[0]);

	for (i=1; i<(M+1)/2-1; i++)
		H[i] = 0.5*h[i*2] + 0.25*(h[i*2-1]+h[i*2+1]);

	return;
}