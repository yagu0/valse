#include <R.h>
#include <Rdefines.h>
#include "selectiontotale.h"

SEXP EMGLLF(
	SEXP phiInit_,
	SEXP rhoInit_,
	SEXP piInit_,
	SEXP gamInit_,
	SEXP mini_,
	SEXP maxi_,
	SEXP gamma_,
	SEXP glambda_,
	SEXP X_,
	SEXP Y_,
	SEXP seuil_,
	SEXP tau_
) {
	// Get matrices dimensions
	SEXP dimX = getAttrib(X_, R_DimSymbol);
	int n = INTEGER(dimX)[0];
	int p = INTEGER(dimX)[1];
	SEXP dimRho = getAttrib(rhoInit_, R_DimSymbol);
	int m = INTEGER(dimRho)[0];
	int k = INTEGER(dimRho)[2];
	int L = INTEGER(getAttrib(glambda_, R_LengthSymbol))[0];

	////////////
	// INPUTS //
	////////////

	// get scalar parameters
	int mini = INTEGER_VALUE(mini_);
	int maxi = INTEGER_VALUE(maxi_);
	double gamma = NUMERIC_VALUE(gamma_);
	double seuil = NUMERIC_VALUE(seuil_);
	double tau = NUMERIC_VALUE(tau_);

	// Get pointers from SEXP arrays ; WARNING: by columns !
	double* piInit = REAL(phiInit_);
	double* rhoInit = REAL(rhoInit_);
	double* piInit = REAL(piInit_);
	double* gamInit = REAL(gamInit_);
	double* glambda = REAL(glambda_);
	double* X = REAL(X_);
	double* Y = REAL(Y_);

	/////////////
	// OUTPUTS //
	/////////////

	int Size = pow(rangmax-rangmin+1,k);
	SEXP A1, A2, rho, pi, dimA;
	PROTECT(dimA = allocVector(INTSXP, 3));
	int* pDimA = INTEGER(dimA);
	pDimA[0] = p; pDimA[1] = m+1; pDimA[2] = L;
	PROTECT(A1 = allocArray(REALSXP, dimA));
	PROTECT(A2 = allocArray(REALSXP, dimA));
	PROTECT(rho = allocArray(REALSXP, dimRho);
	PROTECT(pi = allocMatrix(REALSXP, k, L));
	double *pA1=REAL(A1), *pA2=REAL(A2), *pRho=REAL(rho), *pPi=REAL(pi);

	/////////////////////////////
	// Call to selectiontotale //
	/////////////////////////////

	selectiontotale(
		phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,
		pA1,pA2,pRho,pPi,
		n,p,m,k,L);

	// Build list from OUT params and return it
	SEXP listParams, listNames;
	PROTECT(listParams = allocVector(VECSXP, 4));
	char* lnames[4] = { "A1", "A2", "rho", "pi" }; //lists labels
	PROTECT(listNames = allocVector(STRSXP, 4));
	for (int i=0; i<4; i++)
		SET_STRING_ELT(listNames,i,mkChar(lnames[i]));
	setAttrib(listParams, R_NamesSymbol, listNames);
	SET_ARRAY_ELT(listParams, 0, A1);
	SET_ARRAY_ELT(listParams, 1, A2);
	SET_ARRAY_ELT(listParams, 2, rho);
	SET_MATRIX_ELT(listParams, 3, pi);

	UNPROTECT(7);
	return listParams;
}
