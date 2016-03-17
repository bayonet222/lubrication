#include <stdio.h>
#include <string.h>
#include <math.h>

#define PI 3.141592653

double norm_inf_dif(double *a, double *b, int n);
double norm_inf(double *a, int n);
double norm2_dif(double *a, double *b, int n);

void print_to_file(double *data, int n, int label);
void copyv(double *a, double *b, int n, double alpha);

void main()
{
	const int nx = 1800;
	
	const double dx = 1.0/(double)nx;
	const double dt = 0.2*dx;
	const double k = 2.0*pow(dx,2.0)/dt;

	const	double S = 0.0;
	const double t0 = 0.0, TF = 0.76;

	double *pn, *pk, *pkm, *thetan, *thetak, *thetakm, *cm, *ckm;

	double t;

	double change; const double tol = 5E-8;

	double h, dhdt, h3;

	int ix;

	int PRINT_FREQ = 100000000;

	// MEMORY ALLOCATION
	pk = (double *) malloc(sizeof(double) * nx);
	pkm = (double *) malloc(sizeof(double) * nx);
	thetak = (double *) malloc(sizeof(double) * nx);
	thetakm = (double *) malloc(sizeof(double) * nx);
	cm = (double *) malloc(sizeof(double) * nx);
	ckm = (double *) malloc(sizeof(double) * nx);

	// INITIALIZATION
	for (ix=0; ix<nx; ix++){
		pk[ix] = 0.025;
		thetak[ix] = 1.0;
	}

	copyv(pkm,pk,nx,1.0);
	copyv(thetakm,thetak,nx,1.0);
	copyv(cm,thetak,nx,0.125*cos(4.0*PI*t0) + 0.375);
	copyv(ckm,cm,nx,1.0);

	int niter = 0;

	FILE *fsigma;
	fsigma = fopen("sigma_ea_1800.dat","w");

	double sigma = 0.0;
	int sigmai = 0;


	// TIME ITERATIONS
	while (t <= TF){

		niter += 1; t = t0 + dt*niter;

		h = 0.125*cos(4.0*PI*t) + 0.375;
		h3 = pow(h,3.0);

		change = tol + 1.0;
		while(change > tol) { //GAUSS-SEIDEL

			copyv(ckm,thetak,nx,h);
			
			for (ix=1;ix<nx-1;ix++){ //LOOP OVER NODES
				if (pkm[ix] > 0.0 || thetakm[ix] >= 1.0){

					pk[ix] = 0.5/h3*(S*dx*(ckm[ix-1]-ckm[ix]) - k*(ckm[ix]-cm[ix]) + h3*(pkm[ix+1] + pkm[ix-1]));

					if (pk[ix] >= 0.0) //cavitation check
						thetak[ix] = 1.0;
					else
						pk[ix] = 0.0;
				}

				if (pk[ix] <= 0.0 || thetak[ix] < 1.0){
					thetak[ix] = 1.0/(k+S*dx)/h*(k*cm[ix] + S*dx*ckm[ix-1] + h3*(pkm[ix+1]-2.0*pkm[ix]+pkm[ix-1]));

					if (thetak[ix] < 1.0 ) //cavitation check
						pk[ix] = 0.0;
					else
						thetak[ix] = 1.0;
				}	
			} //END LOOP OVER NODES

			change = norm2_dif(pkm,pk,nx) + norm2_dif(thetakm,thetak,nx);

			copyv(pkm,pk,nx,1.0); copyv(thetakm,thetak,nx,1.0);

		} //END GAUSS-SEIDEL

		copyv(cm,thetak,nx,h);

		printf("tn %1.3e h(t) %1.3f n %d\n",t,h,niter);

		if (niter % PRINT_FREQ == 0)
			print_to_file(pk, nx, niter);

		sigmai = 0;
		sigma = 0.0;

		for (ix=1;ix<nx-1;ix++)
			if (pk[ix]==0.0 && pk[ix+1] > 0.0)
				sigmai = ix;
		sigma = sigmai*dx;

		fprintf(fsigma, "%1.6e %1.6e\n", t, sigma);
			
	} //END while(t <= TF)

	fclose(fsigma);

	free(thetak); free(thetakm);
	free(pk); free(pkm); free(cm); free(ckm);

} 

double norm_inf_dif(double *a, double *b, int n){
	int i;
	double res = 0.0;

	for(i=0;i<n;i++)
		res = fmax(res, fabs(a[i]-b[i]));
	return res;
}

double norm_inf(double *a, int n){
	int i;
	double res = 0.0;

	for(i=0;i<n;i++)
		res = fmax( res, fabs(a[i]) );
	return res;
}

double norm2_dif(double *a, double *b, int n){
	int i;
	double res = 0.0;
	
	for(i=0;i<n;i++)
		res += pow(a[i]-b[i],2.0);
	
	return sqrt(res);
}

void print_to_file(double *data, int n, int label){
	FILE *pfile;
	char filename[20];

	int i;

	snprintf(filename, 20, "res%04d.dat",label);

	pfile = fopen(filename,"w");

	for(i=0;i<n;i++)
		fprintf(pfile, "%1.6e\n", data[i]);

	fclose(pfile);
}

void copyv(double *a, double *b, int n, double alpha){
	int i;

	for(i=0;i<n;i++)
		a[i] = alpha*b[i];
}
