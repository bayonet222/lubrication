#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int N, M, i ,j;
	double *S, *A, *B;
	
	N = mxGetN(prhs[0]);
	M = mxGetM(prhs[0]);
	
	plhs[0] = mxCreateDoubleMatrix(N, M, mxREAL);
	S = mxGetPr(plhs[0]);
	A = mxGetPr(prhs[0]);
	B = mxGetPr(prhs[1]);
	
	for (i=0; i<N; i++)
		for(j=0; j<M; j++)
			S[i+j*N] = A[i+j*N] + B[i+j*N];
	return;
}