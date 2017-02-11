#include <R.h>
#include <Rdefines.h>
#include "constructionModelesLassoRank.h"

SEXP constructionModelesLassoRank(
	SEXP Pi_,
	SEXP Rho_,
	SEXP mini_,
	SEXP maxi_,
	SEXP X_,
	SEXP Y_,
	SEXP tau_,
	SEXP A1_,
	SEXP rangmin_,
	SEXP rangmax_
) {
	// Get matrices dimensions
	SEXP dimX = getAttrib(X_, R_DimSymbol);
	int n = INTEGER(dimX)[0];
	int p = INTEGER(dimX)[1];
	SEXP dimRho = getAttrib(Rho_, R_DimSymbol);
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
	int* A1 = INTEGER(A1_);

	/////////////
	// OUTPUTS //
	/////////////

	int Size = pow(rangmax-rangmin+1,k);
	SEXP phi, llh, dimPhi;
	PROTECT(dimPhi = allocVector(INTSXP, 4));
	int* pDimPhi = INTEGER(dimPhi);
	pDimPhi[0] = p; pDimPhi[1] = m; pDimPhi[2] = k; pDimPhi[3] = L*Size;
	PROTECT(phi = allocArray(REALSXP, dimPhi));
	PROTECT(llh = allocMatrix(REALSXP, L*Size, 2));
	double *pPhi=REAL(phi), *pllh=REAL(llh);

	//////////////////////////////////////////
	// Call to constructionModelesLassoRank //
	//////////////////////////////////////////

	constructionModelesLassoRank_core(
		Pi,Rho,mini,maxi,X,Y,tau,A1,rangmin,rangmax,
		pPhi,pllh,
		n,p,m,k,L);

	// Build list from OUT params and return it
	SEXP listParams, listNames;
	PROTECT(listParams = allocVector(VECSXP, 2));
	char* lnames[2] = {"phi", "llh"}; //lists labels
	PROTECT(listNames = allocVector(STRSXP,2));
	for (int i=0; i<2; i++)
		SET_STRING_ELT(listNames,i,mkChar(lnames[i]));
	setAttrib(listParams, R_NamesSymbol, listNames);
	SET_VECTOR_ELT(listParams, 0, phi);
	SET_VECTOR_ELT(listParams, 1, llh);

	UNPROTECT(5);
	return listParams;
}
