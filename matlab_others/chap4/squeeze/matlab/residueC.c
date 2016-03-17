#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int M, i;
	double *b, *x, *res, dx;
	
	M = mxGetM(prhs[0]);
	
	plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
	res = mxGetPr(plhs[0]);
	
	b = mxGetPr(prhs[0]);
	x = mxGetPr(prhs[1]);
	dx = mxGetScalar(prhs[2]);
	
	for (i=1; i<M-1; i++)
		res[i] = b[i] - 1/pow(dx,2)*(x[i-1] - 2*x[i] + x[i+1]);

	return;
}