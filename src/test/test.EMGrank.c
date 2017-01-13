#include "EMGrank.h"
#include "utils.h"

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
	////////////

	float* Rho = readArray_real("Rho");
	float* Pi = readArray_real("Pi");
	int mini = read_int("mini");
	int maxi = read_int("maxi");
	float* X = readArray_real("X");
	float* Y = readArray_real("Y");
	float tau = read_real("tau");
	int* rank = readArray_int("rank");

	/////////////
	// OUTPUTS //
	/////////////

	float* phi = (float*)malloc(p*m*k*sizeof(float));
	float* LLF = (float*)malloc(1*sizeof(float));

	//////////////////////////
	// Main call to EMGrank //
	//////////////////////////

	EMGrank(Pi,Rho,mini,maxi,X,Y,tau,rank,
		phi,LLF,
		n,p,m,k);

	free(Rho);
	free(Pi);
	free(X);
	free(Y);
	free(rank);

	// Compare to reference outputs
	float* ref_phi = readArray_real("phi");
	compareArray_real("phi", phi, ref_phi, p*m*k);
	free(phi);
	free(ref_phi);

	// LLF
	float* ref_LLF = readArray_real("LLF");
	compareArray_real("LLF", LLF, ref_LLF, 1);
	free(LLF);
	free(ref_LLF);

	return 0;
}
