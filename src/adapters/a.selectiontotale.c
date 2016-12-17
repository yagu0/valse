#include <R.h>
#include <Rdefines.h>
#include "selectiontotale.h"

SEXP selectiontotale(
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
	SEXP dimRhoInit = getAttrib(rhoInit_, R_DimSymbol);
	int m = INTEGER(dimRhoInit)[0];
	int k = INTEGER(dimRhoInit)[2];
	int L = length(glambda_);

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
	double* phiInit = REAL(phiInit_);
	double* rhoInit = REAL(rhoInit_);
	double* piInit = REAL(piInit_);
	double* gamInit = REAL(gamInit_);
	double* glambda = REAL(glambda_);
	double* X = REAL(X_);
	double* Y = REAL(Y_);

	/////////////
	// OUTPUTS //
	/////////////

	SEXP A1, A2, rho, pi, dimA, dimRho;
	PROTECT(dimA = allocVector(INTSXP, 3));
	int* pDimA = INTEGER(dimA);
	pDimA[0] = p; pDimA[1] = m+1; pDimA[2] = L;
	PROTECT(A1 = allocArray(INTSXP, dimA));
	PROTECT(A2 = allocArray(INTSXP, dimA));
	PROTECT(dimRho = allocVector(INTSXP, 4));
	int* pDimRho = INTEGER(dimRho);
	pDimRho[0] = m; pDimRho[1] = m; pDimRho[2] = k; pDimRho[3] = L;
	PROTECT(rho = allocArray(REALSXP, dimRho));
	PROTECT(pi = allocMatrix(REALSXP, k, L));
	int *pA1=INTEGER(A1), *pA2=INTEGER(A2);
	double *pRho=REAL(rho), *pPi=REAL(pi);

	/////////////////////////////
	// Call to selectiontotale //
	/////////////////////////////

	selectiontotale_core(
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
	SET_VECTOR_ELT(listParams, 0, A1);
	SET_VECTOR_ELT(listParams, 1, A2);
	SET_VECTOR_ELT(listParams, 2, rho);
	SET_VECTOR_ELT(listParams, 3, pi);

	UNPROTECT(7);
	return listParams;
}
