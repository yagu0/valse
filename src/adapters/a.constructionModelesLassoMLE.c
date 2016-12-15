#include <R.h>
#include <Rdefines.h>
#include "sources/EMGLLF.h"

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
	SEXP tau_,
	SEXP A1_,
	SEXP A2_
) {
	// Get matrices dimensions
	int n = INTEGER(getAttrib(X_, R_DimSymbol))[0];
	SEXP dim = getAttrib(phiInit_, R_DimSymbol)
	int p = INTEGER(dim)[0];
	int m = INTEGER(dim)[1];
	int k = INTEGER(dim)[2];
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
	double* phiInit = REAL(phiInit_);
	double* rhoInit = REAL(rhoInit_);
	double* piInit = REAL(piInit_);
	double* gamInit = REAL(gamInit_);
	double* glambda = REAL(glambda_);
	double* X = REAL(X_);
	double* Y = REAL(Y_);
	double* A1 = REAL(A1_);
	double* A2 = REAL(A2_);

	/////////////
	// OUTPUTS //
	/////////////

	SEXP phi, rho, pi, lvraisemblance, dimPhi, dimRho;
	PROTECT(dimPhi = allocVector(INTSXP, 4));
	int* pDimPhi = INTEGER(dimPhi);
	pDimPhi[0] = p; pDimPhi[1] = m; pDimPhi[2] = k; pDimPhi[3] = L;
	PROTECT(dimRho = allocVector(INTSXP, 4));
	int* pDimRho = INTEGER(dimRho);
	pDimRho[0] = m; pDimRho[1] = m; pDimRho[2] = k; pDimRho[3] = L;
	PROTECT(phi = allocArray(REALSXP, dimPhi));
	PROTECT(rho = allocArray(REALSXP, dimRho));
	PROTECT(pi = allocMatrix(REALSXP, k, L));
	PROTECT(lvraisemblance = allocMatrix(REALSXP, L, 2));
	double* pPhi=REAL(phi), pRho=REAL(rho), pPi=REAL(pi), pLvraisemblance=REAL(lvraisemblance);

	/////////////////////////////////////////
	// Call to constructionModelesLassoMLE //
	/////////////////////////////////////////

	constructionModelesLassoMLE(
		phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda,X,Y,seuil,tau,A1,A2,
		pPhi,pRho,pPi,pLvraisemblance,
		n,p,m,k,L);

	// Build list from OUT params and return it
	SEXP listParams, listNames;
	PROTECT(listParams = allocVector(VECSXP, 4));
	char* lnames[4] = {"phi", "rho", "pi", "lvraisemblance"}; //lists labels
	PROTECT(listNames = allocVector(STRSXP,4));
	for (int i=0; i<4; i++)
		SET_STRING_ELT(listNames,i,mkChar(lnames[i]));
	setAttrib(listParams, R_NamesSymbol, listNames);
	SET_ARRAY_ELT(listParams, 0, phi);
	SET_ARRAY_ELT(listParams, 1, rho);
	SET_MATRIX_ELT(listParams, 2, pi);
	SET_VECTOR_ELT(listParams, 3, lvraisemblance);

	UNPROTECT(8);
	return listParams;
}
