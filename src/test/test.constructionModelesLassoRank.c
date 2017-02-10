#include "constructionModelesLassoRank.h"
#include "test_utils.h"

int main(int argc, char** argv)
{
	// read dimensions
	const int nbDims = 5;
	int* dimensions = readArray_int("dimensions",&nbDims,1);
	int n = dimensions[0];
	int p = dimensions[1];
	int m = dimensions[2];
	int k = dimensions[3];
	int L = dimensions[4];
	free(dimensions);
	int lengthOne = 1;

	////////////
	// INPUTS //
	////////////

	// piInit
	const int dimPi[] = {k, L};
	float* Pi = readArray_real("Pi",dimPi,2);

	// rhoInit
	const int dimRho[] = {m, m, k, L};
	float* Rho = readArray_real("Rho",dimRho,4);

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
	float* X = readArray_real("X",dimX,2);

	// Y
	const int dimY[] = {n, m};
	float* Y = readArray_real("Y",dimY,2);

	// tau
	float* ptau = readArray_real("tau",&lengthOne,1);
	float tau = *ptau;
	free(ptau);

	// A1
	const int dimA[] = {p, L};
	int* A1 = readArray_int("A1",dimA,2);

	// rangmin
	int* prangmin = readArray_int("rangmin",&lengthOne,1);
	int rangmin = *prangmin;
	free(prangmin);

	// rangmax
	int* prangmax = readArray_int("rangmax",&lengthOne,1);
	int rangmax = *prangmax;
	free(prangmax);
	
	/////////////
	// OUTPUTS //
	/////////////

	// phi
	int Size = (int)pow(rangmax-rangmin+1, k);
	const int dimPhi[] = {p, m, k, L*Size};
	float* phi = (float*)malloc(dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]*sizeof(float));

	// lvraisemblance
	const int dimLvraisemblance[] = {L*Size, 2};
	float* lvraisemblance = (float*)malloc(dimLvraisemblance[0]*dimLvraisemblance[1]*sizeof(float));

	//////////////////////////////////////////////
	// Main call to constructionModelesLassoMLE //
	//////////////////////////////////////////////

	constructionModelesLassoRank(
		Pi,Rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax,
		phi,lvraisemblance,
		n,p,m,k,L);
	
	free(Rho);
	free(Pi);
	free(X);
	free(Y);
	free(A1);
	
	// Compare to reference outputs
	float* ref_phi = readArray_real("phi",dimPhi, 4);
	compareArray_real("phi", phi, ref_phi, dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]);
	free(phi);
	free(ref_phi);
	
	// lvraisemblance
	float* ref_lvraisemblance = readArray_real("lvraisemblance",dimLvraisemblance,2);
	compareArray_real("lvraisemblance", lvraisemblance, ref_lvraisemblance, dimLvraisemblance[0]*dimLvraisemblance[1]);
	free(lvraisemblance);
	free(ref_lvraisemblance);
	
	return 0;
}
