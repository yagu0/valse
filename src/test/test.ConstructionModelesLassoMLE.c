#include "constructionModelesLassoMLE.h"
#include "test_utils.h"

int main(int argc, char** argv)
{
	int* dimensions = readArray_int("dimensions");
	int n = dimensions[0];
	int p = dimensions[1];
	int m = dimensions[2];
	int k = dimensions[3];
	int L = dimensions[4];
	free(dimensions);

	////////////
	// INPUTS //
	Real* phiInit = readArray_real("phiInit");
	Real* rhoInit = readArray_real("rhoInit");
	Real* piInit = readArray_real("piInit");
	Real* gamInit = readArray_real("gamInit");
	int mini = read_int("mini");
	int maxi = read_int("maxi");
	Real gamma = read_real("gamma");
	Real* glambda = readArray_real("glambda");
	Real* X = readArray_real("X");
	Real* Y = readArray_real("Y");
	Real seuil = read_real("seuil");
	Real tau = read_real("tau");
	int* A1 = readArray_int("A1");
	int* A2 = readArray_int("A2");
	////////////

	/////////////
	// OUTPUTS //
	Real* phi = (Real*)malloc(p*m*k*L*sizeof(Real));
	Real* rho = (Real*)malloc(m*m*k*L*sizeof(Real));
	Real* pi = (Real*)malloc(k*L*sizeof(Real));
	Real* llh = (Real*)malloc(L*2*sizeof(Real));
	/////////////

	/////////////////////////////////////////
	// Call to constructionModelesLassoMLE //
	constructionModelesLassoMLE(
		phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,A1,A2,
		phi,rho,pi,llh,
		n,p,m,k,L);
	/////////////////////////////////////////

	free(phiInit);
	free(rhoInit);
	free(piInit);
	free(gamInit);
	free(X);
	free(Y);
	free(A1);
	free(A2);
	free(glambda);

	// Compare to reference outputs
	Real* ref_phi = readArray_real("phi",dimPhi,4);
	compareArray_real("phi", phi, ref_phi, p*m*k*L);
	free(phi);
	free(ref_phi);

	Real* ref_rho = readArray_real("rho",dimRho,4);
	compareArray_real("rho", rho, ref_rho, m*m*k*L);
	free(rho);
	free(ref_rho);

	Real* ref_pi = readArray_real("pi",dimPi,2);
	compareArray_real("pi", pi, ref_pi, k*L);
	free(pi);
	free(ref_pi);

	Real* ref_llh = readArray_real("llh",dimllh,2);
	compareArray_real("llh", llh, ref_llh, L*2);
	free(llh);
	free(ref_llh);

	return 0;
}
