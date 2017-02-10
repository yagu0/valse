#include "selectiontotale.h"
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

	// phiInit
	const int dimPhiInit[] = {p, m, k};
	float* phiInit = readArray_real("phiInit",dimPhiInit,3);

	// rhoInit
	const int dimRhoInit[] = {m, m, k};
	float* rhoInit = readArray_real("rhoInit",dimRhoInit,3);

	// piInit
	float* piInit = readArray_real("piInit",&k,1);

	// gamInit
	const int dimGamInit[] = {n, k};
	float* gamInit = readArray_real("gamInit",dimGamInit, 2);

	// min number of iterations
	int* pmini = readArray_int("mini",&lengthOne,1);
	int mini = *pmini;
	free(pmini);
	
	// max number of iterations
	int* pmaxi = readArray_int("maxi",&lengthOne,1);
	int maxi = *pmaxi;
	free(pmaxi);
	
	// gamma
	float* pgamma = readArray_real("gamma",&lengthOne,1);
	float gamma = *pgamma;
	free(pgamma);
	
	// lambda
	float* glambda = readArray_real("glambda",&L,1);

	// X
	const int dimX[] = {n, p};
	float* X = readArray_real("X",dimX,2);

	// Y
	const int dimY[] = {n, m};
	float* Y = readArray_real("Y",dimY,2);

	// seuil
	float* pseuil = readArray_real("seuil",&lengthOne,1);
	float seuil = *pseuil;
	free(pseuil);
	
	// tau
	float* ptau = readArray_real("tau",&lengthOne,1);
	float tau = *ptau;
	free(ptau);
	
	/////////////
	// OUTPUTS //
	/////////////

	// A1
	const int dimA[] = {p, m+1, L};
	int* A1 = (int*)malloc(dimA[0]*dimA[1]*dimA[2]*sizeof(int));

	// A2
	int* A2 = (int*)malloc(dimA[0]*dimA[1]*dimA[2]*sizeof(int));

	// rho
	const int dimRho[] = {m, m, k, L};
	float* Rho = (float*)malloc(dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]*sizeof(float));

	// pi
	const int dimPi[] = {k, L};
	float* Pi = (float*)malloc(dimPi[0]*dimPi[1]*sizeof(float));

	//////////////////////////////////////////////
	// Main call to constructionModelesLassoMLE //
	//////////////////////////////////////////////

	selectiontotale(
		phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,
		A1,A2,Rho,Pi,
		n,p,m,k,L);
	
	// free input pointers
	free(phiInit);
	free(rhoInit);
	free(piInit);
	free(gamInit);
	free(glambda);
	free(X);
	free(Y);
	
	// Compare to reference outputs
	int* ref_A1 = readArray_int("A1",dimA, 3);
	compareArray_int("A1", A1, ref_A1, dimA[0]*dimA[1]*dimA[2]);
	free(A1);
	free(ref_A1);
	
	// A2
	int* ref_A2 = readArray_int("A2",dimA, 3);
	compareArray_int("A2", A2, ref_A2, dimA[0]*dimA[1]*dimA[2]);
	free(A2);
	free(ref_A2);
	
	// Rho
	float* ref_Rho = readArray_real("Rho",dimRho, 4);
	compareArray_real("Rho", Rho, ref_Rho, dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]);
	free(Rho);
	free(ref_Rho);
	
	// Pi
	float* ref_Pi = readArray_real("Pi",dimPi, 2);
	compareArray_real("Pi", Pi, ref_Pi, dimPi[0]*dimPi[1]);
	free(Pi);
	free(ref_Pi);
	
	return 0;
}
