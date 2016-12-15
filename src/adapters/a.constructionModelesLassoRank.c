#include <R.h>
#include <Rdefines.h>
#include "sources/EMGLLF.h"

SEXP EMGLLF(
	SEXP Pi_,
	SEXP Rho_,
	SEXP mini_,
	SEXP maxi_,
	SEXP X_,
	SEXP Y_,
	SEXP tau_,
	SEXP A1_,
	SEXP rangmin_,
	SEXP rangmax
) {
	// Get matrices dimensions
	SEXP dimX = getAttrib(X_, R_DimSymbol);
	int n = INTEGER(dimX)[0];
	int p = INTEGER(dimX)[1];
	SEXP dimRho = getAttrib(Rho_, R_DimSymbol)
	int m = INTEGER(dimRho)[0];
	int k = INTEGER(dimRho)[2];
	int L = INTEGER(getAttrib(A1_, R_DimSymbol))[1];

	////////////
	// INPUTS //
	////////////

	// get scalar parameters
	int mini = INTEGER_VALUE(mini_);
	int maxi = INTEGER_VALUE(maxi_);
	double tau = NUMERIC_VALUE(tau_);
	double rangmin = NUMERIC_VALUE(rangmin_);
	double rangmax = NUMERIC_VALUE(rangmax_);

	// Get pointers from SEXP arrays ; WARNING: by columns !
	double* Pi = REAL(Pi_);
	double* Rho = REAL(Rho_);
	double* X = REAL(X_);
	double* Y = REAL(Y_);
	double* A1 = REAL(A1_);

	/////////////
	// OUTPUTS //
	/////////////

	int Size = pow(rangmax-rangmin+1,k);
	SEXP phi, lvraisemblance, dimPhi;
	PROTECT(dimPhi = allocVector(INTSXP, 4));
	int* pDimPhi = INTEGER(dimPhi);
	pDimPhi[0] = p; pDimPhi[1] = m; pDimPhi[2] = k; pDimPhi[3] = L*Size;
	PROTECT(phi = allocArray(REALSXP, dimPhi));
	PROTECT(lvraisemblance = allocMatrix(REALSXP, L*Size, 2));
	double* pPhi=REAL(phi), pLvraisemblance=REAL(lvraisemblance);

	//////////////////////////////////////////
	// Call to constructionModelesLassoRank //
	//////////////////////////////////////////

	constructionModelesLassoRank(
		Pi,Rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax,
		pPhi,pLvraisemblance,
		n,p,m,k,L);

	// Build list from OUT params and return it
	SEXP listParams, listNames;
	PROTECT(listParams = allocVector(VECSXP, 2));
	char* lnames[2] = {"phi", "lvraisemblance"}; //lists labels
	PROTECT(listNames = allocVector(STRSXP,2));
	for (int i=0; i<2; i++)
		SET_STRING_ELT(listNames,i,mkChar(lnames[i]));
	setAttrib(listParams, R_NamesSymbol, listNames);
	SET_ARRAY_ELT(listParams, 0, phi);
	SET_VECTOR_ELT(listParams, 1, lvraisemblance);

	UNPROTECT(5);
	return listParams;
}
