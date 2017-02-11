#include "selectiontotale.h"
#include "test_utils.h"
#include <stdlib.h>

int main(int argc, char** argv)
{
	int* dimensions = readArray_int("dimensions");
	int n = dimensions[0];
	int p = dimensions[1];
	int m = dimensions[2];
	int k = dimensions[3];
	int L = dimensions[4];
	free(dimensions);

	////////////
	// INPUTS //
	Real* phiInit = readArray_real("phiInit");
	Real* rhoInit = readArray_real("rhoInit");
	Real* piInit = readArray_real("piInit");
	Real* gamInit = readArray_real("gamInit");
	int mini = read_int("mini");
	int maxi = read_int("maxi");
	Real gamma = read_real("gamma");
	Real* glambda = readArray_real("glambda");
	Real* X = readArray_real("X");
	Real* Y = readArray_real("Y");
	Real seuil = read_real("seuil");
	Real tau = read_real("tau");
	////////////

	/////////////
	// OUTPUTS //
	int* A1 = (int*)malloc(p*(m+1)*L*sizeof(int));
	int* A2 = (int*)malloc(p*(m+1)*L*sizeof(int));
	Real* Rho = (Real*)malloc(m*m*k*L*sizeof(Real));
	Real* Pi = (Real*)malloc(k*L*sizeof(Real));
	/////////////

	/////////////////////////////////////////
	// Call to constructionModelesLassoMLE //
	selectiontotale_core(
		phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,
		A1,A2,Rho,Pi,
		n,p,m,k,L);
	/////////////////////////////////////////

	free(phiInit);
	free(rhoInit);
	free(piInit);
	free(gamInit);
	free(glambda);
	free(X);
	free(Y);

	// Compare to reference outputs
	int* ref_A1 = readArray_int("A1");
	compareArray_int("A1", A1, ref_A1, p*(m+1)*L);
	free(A1);
	free(ref_A1);

	int* ref_A2 = readArray_int("A2");
	compareArray_int("A2", A2, ref_A2, p*(m+1)*L);
	free(A2);
	free(ref_A2);

	Real* ref_Rho = readArray_real("Rho");
	compareArray_real("Rho", Rho, ref_Rho, m*m*k*L);
	free(Rho);
	free(ref_Rho);

	Real* ref_Pi = readArray_real("Pi");
	compareArray_real("Pi", Pi, ref_Pi, k*L);
	free(Pi);
	free(ref_Pi);

	return 0;
}
