#include "constructionModelesLassoRank.h"
#include "test_utils.h"
#include <stdlib.h>
#include <math.h>

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
	Real* pi = readArray_real("pi");
	Real* rho = readArray_real("rho");
	int mini = read_int("mini");
	int maxi = read_int("maxi");
	Real* X = readArray_real("X");
	Real* Y = readArray_real("Y");
	Real tau = read_real("tau");
	int* A1 = readArray_int("A1");
	int rangmin = read_int("rangmin");
	int rangmax = read_int("rangmax");
	////////////

	/////////////
	// OUTPUTS //
	int Size = (int)pow(rangmax-rangmin+1, k);
	Real* phi = (Real*)malloc(p*m*k*L*Size*sizeof(Real));
	Real* llh = (Real*)malloc(L*Size*2*sizeof(Real));
	/////////////

	/////////////////////////////////////////
	// Call to constructionModelesLassoMLE //
	constructionModelesLassoRank_core(
		pi,rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax,
		phi,llh,
		n,p,m,k,L);
	/////////////////////////////////////////

	free(rho);
	free(pi);
	free(X);
	free(Y);
	free(A1);

	// Compare to reference outputs
	Real* ref_phi = readArray_real("phi");
	compareArray_real("phi", phi, ref_phi, p*m*k*L*Size);
	free(phi);
	free(ref_phi);

	Real* ref_llh = readArray_real("llh");
	compareArray_real("llh", llh, ref_llh, L*Size*2);
	free(llh);
	free(ref_llh);

	return 0;
}
