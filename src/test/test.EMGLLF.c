#include "EMGLLF.h"
#include "ioutils.h"

int main(int argc, char** argv)
{
	////////////
	// INPUTS //
	////////////

	Int* dimensions = readArray_int("dimensions");
	mwSize n = dimensions[0];
	mwSize p = dimensions[1];
	mwSize m = dimensions[2];
	mwSize k = dimensions[3];
	free(dimensions);

	Real* phiInit = readArray_real("phiInit");
	Real* rhoInit = readArray_real("rhoInit");
	Real* piInit = readArray_real("piInit");
	Real* gamInit = readArray_real("gamInit");
	Int mini = read_int("mini");
	Int maxi = read_int("maxi");
	Real gamma = read_real("gamma");
	Real lambda = read_real("lambda");
	Real* X = readArray_real("X");
	Real* Y = readArray_real("Y");
	Real tau = read_real("tau");

	/////////////
	// OUTPUTS //
	/////////////

	Real* phi = (Real*)malloc(p*m*k*sizeof(Real));
	Real* rho = (Real*)malloc(m*m*k*sizeof(Real));
	Real* pi = (Real*)malloc(k*sizeof(Real));
	Real* LLF = (Real*)malloc(maxi*sizeof(Real));
	Real* S = (Real*)malloc(p*m*k*sizeof(Real));

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
	Real* ref_phi = readArray_real("phi");
	compareArray_real("phi", phi, ref_phi, p*m*k);
	free(phi);
	free(ref_phi);

	Real* ref_rho = readArray_real("rho");
	compareArray_real("rho", rho, ref_rho, m*m*k);
	free(rho);
	free(ref_rho);

	Real* ref_pi = readArray_real("pi");
	compareArray_real("pi", pi, ref_pi, k);
	free(pi);
	free(ref_pi);

	Real* ref_LLF = readArray_real("LLF", maxi);
	compareArray_real("LLF", LLF, ref_LLF);
	free(LLF);
	free(ref_LLF);

	Real* ref_S = readArray_real("S");
	compareArray_real("S", S, ref_S, p*m*k);
	free(S);
	free(ref_S);

	return 0;
}
