#include "EMGrank.h"
#include "ioutils.h"

int main(int argc, char** argv)
{
	// read dimensions
	const int nbDims = 4;
	int* dimensions = readArray_int("dimensions",&nbDims,1);
	int n = dimensions[0];
	int p = dimensions[1];
	int m = dimensions[2];
	int k = dimensions[3];
	free(dimensions);
	int lengthOne = 1;

	////////////
	// INPUTS //
	////////////

	// Rho
	const int dimRho[] = {m, m, k};
	float* Rho = readArray_real("Rho",dimRho,3);

	// Pi
	float* Pi = readArray_real("Pi",&k,1);

	// min number of iterations
	int* pmini = readArray_int("mini",&lengthOne,1);
	int mini = *pmini;
	free(pmini);

	// max number of iterations
	int* pmaxi = readArray_int("maxi",&lengthOne,1);
	int maxi = *pmaxi;
	free(pmaxi);

	// X
	const int dimX[] = {n, p};
	float* X = readArray_real("X",dimX, 2);

	// Y
	const int dimY[] = {n, m};
	float* Y = readArray_real("Y",dimY, 2);

	// tau
	float* ptau = readArray_real("tau",&lengthOne,1);
	float tau = *ptau;
	free(ptau);

	// tau
	int* rank = readArray_int("rank",&k,1);

	/////////////
	// OUTPUTS //
	/////////////

	// phi
	const int dimPhi[] = {p, m, k};
	float* phi = (float*)malloc(dimPhi[0]*dimPhi[1]*dimPhi[2]*sizeof(float));

	// LLF
	float* LLF = (float*)malloc(1*sizeof(float));

	//////////////////////////
	// Main call to EMGrank //
	//////////////////////////

	EMGrank(Pi,Rho,mini,maxi,X,Y,tau,rank,
		phi,LLF,
		n,p,m,k);
	
	// free input pointers
	free(Rho);
	free(Pi);
	free(X);
	free(Y);
	free(rank);
	
	// Compare to reference outputs
	float* ref_phi = readArray_real("phi",dimPhi, 3);
	compareArray_real("phi", phi, ref_phi, dimPhi[0]*dimPhi[1]*dimPhi[2]);
	free(phi);
	free(ref_phi);
	
	// LLF
	float* ref_LLF = readArray_real("LLF",&lengthOne,1);
	compareArray_real("LLF", LLF, ref_LLF, 1);
	free(LLF);
	free(ref_LLF);
	
	return 0;
}
