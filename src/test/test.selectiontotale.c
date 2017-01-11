#include "selectiontotale.h"
#include "ioutils.h"

mwSize main(mwSize argc, char** argv)
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
	Real* gamInit = readArray_real("gamInit",dimGamInit, 2);

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
	
	/////////////
	// OUTPUTS //
	/////////////

	// A1
	const mwSize dimA[] = {p, m+1, L};
	Int* A1 = (Int*)malloc(dimA[0]*dimA[1]*dimA[2]*sizeof(Int));

	// A2
	Int* A2 = (Int*)malloc(dimA[0]*dimA[1]*dimA[2]*sizeof(Int));

	// rho
	const mwSize dimRho[] = {m, m, k, L};
	Real* Rho = (Real*)malloc(dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]*sizeof(Real));

	// pi
	const mwSize dimPi[] = {k, L};
	Real* Pi = (Real*)malloc(dimPi[0]*dimPi[1]*sizeof(Real));

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
	Int* ref_A1 = readArray_int("A1",dimA, 3);
	compareArray_int("A1", A1, ref_A1, dimA[0]*dimA[1]*dimA[2]);
	free(A1);
	free(ref_A1);
	
	// A2
	Int* ref_A2 = readArray_int("A2",dimA, 3);
	compareArray_int("A2", A2, ref_A2, dimA[0]*dimA[1]*dimA[2]);
	free(A2);
	free(ref_A2);
	
	// Rho
	Real* ref_Rho = readArray_real("Rho",dimRho, 4);
	compareArray_real("Rho", Rho, ref_Rho, dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]);
	free(Rho);
	free(ref_Rho);
	
	// Pi
	Real* ref_Pi = readArray_real("Pi",dimPi, 2);
	compareArray_real("Pi", Pi, ref_Pi, dimPi[0]*dimPi[1]);
	free(Pi);
	free(ref_Pi);
	
	return 0;
}
