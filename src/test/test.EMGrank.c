#include "EMGrank.h"
#include "ioutils.h"

int main(int argc, char** argv)
{
	// read dimensions
	const Int nbDims = 4;
	Int* dimensions = readArray_int("dimensions",&nbDims,1);
	mwSize n = dimensions[0];
	mwSize p = dimensions[1];
	mwSize m = dimensions[2];
	mwSize k = dimensions[3];
	free(dimensions);
	mwSize lengthOne = 1;

	////////////
	// INPUTS //
	////////////

	// Rho
	const mwSize dimRho[] = {m, m, k};
	Real* Rho = readArray_real("Rho",dimRho,3);

	// Pi
	Real* Pi = readArray_real("Pi",&k,1);

	// min number of iterations
	Int* pmini = readArray_int("mini",&lengthOne,1);
	Int mini = *pmini;
	free(pmini);

	// max number of iterations
	Int* pmaxi = readArray_int("maxi",&lengthOne,1);
	Int maxi = *pmaxi;
	free(pmaxi);

	// X
	const mwSize dimX[] = {n, p};
	Real* X = readArray_real("X",dimX, 2);

	// Y
	const mwSize dimY[] = {n, m};
	Real* Y = readArray_real("Y",dimY, 2);

	// tau
	Real* ptau = readArray_real("tau",&lengthOne,1);
	Real tau = *ptau;
	free(ptau);

	// tau
	Int* rank = readArray_int("rank",&k,1);

	/////////////
	// OUTPUTS //
	/////////////

	// phi
	const mwSize dimPhi[] = {p, m, k};
	Real* phi = (Real*)malloc(dimPhi[0]*dimPhi[1]*dimPhi[2]*sizeof(Real));

	// LLF
	Real* LLF = (Real*)malloc(1*sizeof(Real));

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
	Real* ref_phi = readArray_real("phi",dimPhi, 3);
	compareArray_real("phi", phi, ref_phi, dimPhi[0]*dimPhi[1]*dimPhi[2]);
	free(phi);
	free(ref_phi);
	
	// LLF
	Real* ref_LLF = readArray_real("LLF",&lengthOne,1);
	compareArray_real("LLF", LLF, ref_LLF, 1);
	free(LLF);
	free(ref_LLF);
	
	return 0;
}
