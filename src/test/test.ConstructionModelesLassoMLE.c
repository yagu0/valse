#include "constructionModelesLassoMLE.h"
#include "ioutils.h"

int main(int argc, char** argv)
{
	// read dimensions
	const Int nbDims = 5;
	Int* dimensions = readArray_int("dimensions",&nbDims,1);
	mwSize n = dimensions[0];
	mwSize p = dimensions[1];
	mwSize m = dimensions[2];
	mwSize k = dimensions[3];
	mwSize L = dimensions[4];
	free(dimensions);
	mwSize lengthOne = 1;

	////////////
	// INPUTS //
	////////////

	// phiInit
	const mwSize dimPhiInit[] = {p, m, k};
	Real* phiInit = readArray_real("phiInit",dimPhiInit,3);

	// rhoInit
	const mwSize dimRhoInit[] = {m, m, k};
	Real* rhoInit = readArray_real("rhoInit",dimRhoInit,3);

	// piInit
	Real* piInit = readArray_real("piInit",&k,1);

	// gamInit
	const mwSize dimGamInit[] = {n, k};
	Real* gamInit = readArray_real("gamInit",dimGamInit,2);

	// min number of iterations
	Int* pmini = readArray_int("mini",&lengthOne,1);
	Int mini = *pmini;
	free(pmini);
	
	// max number of iterations
	Int* pmaxi = readArray_int("maxi",&lengthOne,1);
	Int maxi = *pmaxi;
	free(pmaxi);

	// gamma
	Real* pgamma = readArray_real("gamma",&lengthOne,1);
	Real gamma = *pgamma;
	free(pgamma);

	// lambda
	Real* glambda = readArray_real("glambda",&L,1);

	// X
	const mwSize dimX[] = {n, p};
	Real* X = readArray_real("X",dimX,2);

	// Y
	const mwSize dimY[] = {n, m};
	Real* Y = readArray_real("Y",dimY,2);

	// seuil
	Real* pseuil = readArray_real("seuil",&lengthOne,1);
	Real seuil = *pseuil;
	free(pseuil);

	// tau
	Real* ptau = readArray_real("tau",&lengthOne,1);
	Real tau = *ptau;
	free(ptau);
	
	// A1
	const mwSize dimA[] = {p, m+1, L};
	Int* A1 = readArray_int("A1",dimA,3);
	
	// A2
	Int* A2 = readArray_int("A2",dimA,3);
	
	/////////////
	// OUTPUTS //
	/////////////

	// phi
	const mwSize dimPhi[] = {dimPhiInit[0], dimPhiInit[1], dimPhiInit[2], L};
	Real* phi = (Real*)malloc(dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]*sizeof(Real));

	// rho
	const mwSize dimRho[] = {dimRhoInit[0], dimRhoInit[1], dimRhoInit[2], L};
	Real* rho = (Real*)malloc(dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]*sizeof(Real));

	// pi
	const mwSize dimPi[] = {k, L};
	Real* pi = (Real*)malloc(dimPi[0]*dimPi[1]*sizeof(Real));

	// lvraisemblance
	const mwSize dimLvraisemblance[] = {L, 2};
	Real* lvraisemblance = (Real*)malloc(dimLvraisemblance[0]*dimLvraisemblance[1]*sizeof(Real));

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
	Real* ref_phi = readArray_real("phi",dimPhi,4);
	compareArray_real("phi", phi, ref_phi, dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]);
	free(phi);
	free(ref_phi);
	
	// rho
	Real* ref_rho = readArray_real("rho",dimRho,4);
	compareArray_real("rho", rho, ref_rho, dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]);
	free(rho);
	free(ref_rho);
	
	// pi
	Real* ref_pi = readArray_real("pi",dimPi,2);
	compareArray_real("pi", pi, ref_pi, dimPi[0]*dimPi[1]);
	free(pi);
	free(ref_pi);
	
	// lvraisemblance
	Real* ref_lvraisemblance = readArray_real("lvraisemblance",dimLvraisemblance,2);
	compareArray_real("lvraisemblance", lvraisemblance, ref_lvraisemblance, dimLvraisemblance[0]*dimLvraisemblance[1]);
	free(lvraisemblance);
	free(ref_lvraisemblance);

	return 0;
}
