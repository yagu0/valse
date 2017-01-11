#include "constructionModelesLassoRank.h"
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

	// piInit
	const mwSize dimPi[] = {k, L};
	Real* Pi = readArray_real("Pi",dimPi,2);

	// rhoInit
	const mwSize dimRho[] = {m, m, k, L};
	Real* Rho = readArray_real("Rho",dimRho,4);

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
	Real* X = readArray_real("X",dimX,2);

	// Y
	const mwSize dimY[] = {n, m};
	Real* Y = readArray_real("Y",dimY,2);

	// tau
	Real* ptau = readArray_real("tau",&lengthOne,1);
	Real tau = *ptau;
	free(ptau);

	// A1
	const mwSize dimA[] = {p, L};
	Int* A1 = readArray_int("A1",dimA,2);

	// rangmin
	Int* prangmin = readArray_int("rangmin",&lengthOne,1);
	Int rangmin = *prangmin;
	free(prangmin);

	// rangmax
	Int* prangmax = readArray_int("rangmax",&lengthOne,1);
	Int rangmax = *prangmax;
	free(prangmax);
	
	/////////////
	// OUTPUTS //
	/////////////

	// phi
	mwSize Size = (mwSize)pow(rangmax-rangmin+1, k);
	const mwSize dimPhi[] = {p, m, k, L*Size};
	Real* phi = (Real*)malloc(dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]*sizeof(Real));

	// lvraisemblance
	const mwSize dimLvraisemblance[] = {L*Size, 2};
	Real* lvraisemblance = (Real*)malloc(dimLvraisemblance[0]*dimLvraisemblance[1]*sizeof(Real));

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
	Real* ref_phi = readArray_real("phi",dimPhi, 4);
	compareArray_real("phi", phi, ref_phi, dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]);
	free(phi);
	free(ref_phi);
	
	// lvraisemblance
	Real* ref_lvraisemblance = readArray_real("lvraisemblance",dimLvraisemblance,2);
	compareArray_real("lvraisemblance", lvraisemblance, ref_lvraisemblance, dimLvraisemblance[0]*dimLvraisemblance[1]);
	free(lvraisemblance);
	free(ref_lvraisemblance);
	
	return 0;
}
