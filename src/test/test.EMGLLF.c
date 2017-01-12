#include "EMGLLF.h"
#include "utils.h"
#include <stdlib.h>

int main(int argc, char** argv)
{
	////////////
	// INPUTS //
	////////////

	int* dimensions = readArray_int("dimensions");
	int n = dimensions[0];
	int p = dimensions[1];
	int m = dimensions[2];
	int k = dimensions[3];
	free(dimensions);

	float* phiInit = readArray_real("phiInit");
	float* rhoInit = readArray_real("rhoInit");
	float* piInit = readArray_real("piInit");
	float* gamInit = readArray_real("gamInit");
	int mini = read_int("mini");
	int maxi = read_int("maxi");
	float gamma = read_real("gamma");
	float lambda = read_real("lambda");
	float* X = readArray_real("X");
	float* Y = readArray_real("Y");
	float tau = read_real("tau");

	/////////////
	// OUTPUTS //
	/////////////

	float* phi = (float*)malloc(p*m*k*sizeof(float));
	float* rho = (float*)malloc(m*m*k*sizeof(float));
	float* pi = (float*)malloc(k*sizeof(float));
	float* LLF = (float*)malloc(maxi*sizeof(float));
	float* S = (float*)malloc(p*m*k*sizeof(float));

	////////////////////
	// Call to EMGLLF //
	////////////////////

	EMGLLF_core(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,lambda,X,Y,tau,
		phi,rho,pi,LLF,S,
		n,p,m,k);

	free(phiInit);
	free(rhoInit);
	free(piInit);
	free(gamInit);
	free(X);
	free(Y);

	// Compare to reference outputs
	float* ref_phi = readArray_real("phi");
	compareArray_real("phi", phi, ref_phi, p*m*k);
	free(phi);
	free(ref_phi);

	float* ref_rho = readArray_real("rho");
	compareArray_real("rho", rho, ref_rho, m*m*k);
	free(rho);
	free(ref_rho);

	float* ref_pi = readArray_real("pi");
	compareArray_real("pi", pi, ref_pi, k);
	free(pi);
	free(ref_pi);

	float* ref_LLF = readArray_real("LLF", maxi);
	compareArray_real("LLF", LLF, ref_LLF);
	free(LLF);
	free(ref_LLF);

	float* ref_S = readArray_real("S");
	compareArray_real("S", S, ref_S, p*m*k);
	free(S);
	free(ref_S);

	return 0;
}
