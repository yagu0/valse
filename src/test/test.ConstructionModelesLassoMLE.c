#include "constructionModelesLassoMLE.h"
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
	float* gamInit = readArray_real("gamInit",dimGamInit,2);

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
	
	// A1
	const int dimA[] = {p, m+1, L};
	int* A1 = readArray_int("A1",dimA,3);
	
	// A2
	int* A2 = readArray_int("A2",dimA,3);
	
	/////////////
	// OUTPUTS //
	/////////////

	// phi
	const int dimPhi[] = {dimPhiInit[0], dimPhiInit[1], dimPhiInit[2], L};
	float* phi = (float*)malloc(dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]*sizeof(float));

	// rho
	const int dimRho[] = {dimRhoInit[0], dimRhoInit[1], dimRhoInit[2], L};
	float* rho = (float*)malloc(dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]*sizeof(float));

	// pi
	const int dimPi[] = {k, L};
	float* pi = (float*)malloc(dimPi[0]*dimPi[1]*sizeof(float));

	// lvraisemblance
	const int dimLvraisemblance[] = {L, 2};
	float* lvraisemblance = (float*)malloc(dimLvraisemblance[0]*dimLvraisemblance[1]*sizeof(float));

	/////////////////////////////////////////
	// Call to constructionModelesLassoMLE //
	/////////////////////////////////////////

	constructionModelesLassoMLE(
		phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,A1,A2,
		phi,rho,pi,lvraisemblance,
		n,p,m,k,L);
	
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
	float* ref_phi = readArray_real("phi",dimPhi,4);
	compareArray_real("phi", phi, ref_phi, dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]);
	free(phi);
	free(ref_phi);
	
	// rho
	float* ref_rho = readArray_real("rho",dimRho,4);
	compareArray_real("rho", rho, ref_rho, dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]);
	free(rho);
	free(ref_rho);
	
	// pi
	float* ref_pi = readArray_real("pi",dimPi,2);
	compareArray_real("pi", pi, ref_pi, dimPi[0]*dimPi[1]);
	free(pi);
	free(ref_pi);
	
	// lvraisemblance
	float* ref_lvraisemblance = readArray_real("lvraisemblance",dimLvraisemblance,2);
	compareArray_real("lvraisemblance", lvraisemblance, ref_lvraisemblance, dimLvraisemblance[0]*dimLvraisemblance[1]);
	free(lvraisemblance);
	free(ref_lvraisemblance);

	return 0;
}
