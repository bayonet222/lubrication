#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int M, i, iiter, niter, flag;
	double *b, *res, *x0, dx;
	
	M = mxGetM(prhs[0]);
	
	plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
	res = mxGetPr(plhs[0]);
	
	b = mxGetPr(prhs[0]);
	dx = mxGetScalar(prhs[1]);
	
	niter = floor( mxGetScalar(prhs[2]) );
	
	x0 = mxGetPr(prhs[3]);
	flag = floor( mxGetScalar(prhs[4]) );
	
	if (flag!=0)
		for (i=1; i<M-1; i++)
			res[i] = x0[i];
	
	res[0] = x0[0]; res[M-1] = x0[M-1];
	
	for (iiter=0; iiter < niter; iiter++)
		for (i=1; i<M-1; i++)
			res[i] = -0.5*(pow(dx,2)*b[i]-res[i-1]-res[i+1]);
		
	return;
}