#include "EMGrank.h"
#include "test_utils.h"
#include <stdlib.h>

int main(int argc, char** argv)
{
	int* dimensions = readArray_int("dimensions");
	int n = dimensions[0];
	int p = dimensions[1];
	int m = dimensions[2];
	int k = dimensions[3];
	free(dimensions);

	////////////
	// INPUTS //
	Real* Rho = readArray_real("Rho");
	Real* Pi = readArray_real("Pi");
	int mini = read_int("mini");
	int maxi = read_int("maxi");
	Real* X = readArray_real("X");
	Real* Y = readArray_real("Y");
	Real eps = read_real("eps");
	int* rank = readArray_int("rank");
	////////////

	/////////////
	// OUTPUTS //
	Real* phi = (Real*)malloc(p*m*k*sizeof(Real));
	Real* LLF = (Real*)malloc(1*sizeof(Real));
	/////////////

	//////////////////////////
	// Main call to EMGrank //
	EMGrank_core(Pi,Rho,mini,maxi,X,Y,eps,rank,
		phi,LLF,
		n,p,m,k);
	//////////////////////////

	free(Rho);
	free(Pi);
	free(X);
	free(Y);
	free(rank);

	// Compare to reference outputs
	Real* ref_phi = readArray_real("phi");
	compareArray_real("phi", phi, ref_phi, p*m*k);
	free(phi);
	free(ref_phi);

	// LLF
	Real* ref_LLF = readArray_real("LLF");
	compareArray_real("LLF", LLF, ref_LLF, 1);
	free(LLF);
	free(ref_LLF);

	return 0;
}
