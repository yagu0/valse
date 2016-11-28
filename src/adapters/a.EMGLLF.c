#include <R.h>
#include <Rdefines.h>
#include "sources/EMGLLF.h"
#include "sources/utils/io.h"

SEXP EMGLLF(
	SEXP M_, 
	SEXP NIix_, 
	SEXP alpha_, 
	SEXP h_, 
	SEXP epsilon_, 
	SEXP maxiter_, 
	SEXP symmNeighbs_, 
	SEXP trace_
) {
	// get parameters
	double alpha = NUMERIC_VALUE(alpha_);
	double h = NUMERIC_VALUE(h_);
	double epsilon = NUMERIC_VALUE(epsilon_);
	int maxiter = INTEGER_VALUE(maxiter_);
	int symmNeighbs = LOGICAL_VALUE(symmNeighbs_);
	int trace = LOGICAL_VALUE(trace_);

	// extract infos from M and get associate pointer
	SEXP dim = getAttrib(M_, R_DimSymbol);
	int nrow = INTEGER(dim)[0];
	int ncol = INTEGER(dim)[1];
	// M is always given by columns: easier to process in rows
	double* pM = transpose(REAL(M_), nrow, ncol);

	// extract NIix list vectors in a jagged array
	int* lengthNIix = (int*)malloc(nrow*sizeof(int));
	int** NIix = (int**)malloc(nrow*sizeof(int*));
	for (int i=0; i<nrow; i++)
	{
		lengthNIix[i] = LENGTH(VECTOR_ELT(NIix_,i));
		SEXP tmp;
		PROTECT(tmp = AS_INTEGER(VECTOR_ELT(NIix_,i)));
		NIix[i] = (int*)malloc(lengthNIix[i]*sizeof(int));
		for (int j=0; j<lengthNIix[i]; j++)
			NIix[i][j] = INTEGER(tmp)[j];
		UNPROTECT(1);
		// WARNING: R indices start at 1,
		// so we must lower every index right now to avoid future bug
		for (int j=0; j<lengthNIix[i]; j++)
			NIix[i][j]--;
	}

	// Main call to core algorithm
	Parameters params = getVarsWithConvexOptim_core(
		pM, lengthNIix, NIix, nrow, ncol, alpha, h, epsilon, maxiter, symmNeighbs, trace);

	// free neighborhoods parameters arrays
	free(lengthNIix);
	for (int i=0; i<nrow; i++)
		free(NIix[i]);
	free(NIix);

	// copy matrix F into pF for output to R (1D matrices)
	SEXP f;
	PROTECT(f = allocMatrix(REALSXP, nrow, ncol));
	double* pF = REAL(f);
	for (int i=0; i<nrow; i++)
	{
		for (int j=0; j<ncol; j++)
			pF[i+nrow*j] = params.f[i][j];
	}
	// copy theta into pTheta for output to R
	SEXP theta;
	PROTECT(theta = allocVector(REALSXP, nrow));
	double* pTheta = REAL(theta);
	for (int i=0; i<nrow; i++)
		pTheta[i] = params.theta[i];

	// free params.f and params.theta
	free(params.theta);
	for (int i=0; i<nrow; i++)
		free(params.f[i]);
	free(params.f);

	// build return list with f and theta
	SEXP listParams, listNames;
	PROTECT(listParams = allocVector(VECSXP, 2));
	char* lnames[2] = {"f", "theta"}; //lists labels
	PROTECT(listNames = allocVector(STRSXP,2));
	for (int i=0; i<2; i++)
		SET_STRING_ELT(listNames,i,mkChar(lnames[i]));
	setAttrib(listParams, R_NamesSymbol, listNames);
	SET_VECTOR_ELT(listParams, 0, f);
	SET_VECTOR_ELT(listParams, 1, theta);

	UNPROTECT(4);
	return listParams;








	// Get matrices dimensions
	const mwSize n = mxGetDimensions(prhs[8])[0];
	const mwSize p = mxGetDimensions(prhs[0])[0];
	const mwSize m = mxGetDimensions(prhs[0])[1];
	const mwSize k = mxGetDimensions(prhs[0])[2];

	////////////
	// INPUTS //
	////////////

	// phiInit
	const mwSize* dimPhiInit = mxGetDimensions(prhs[0]);
	Real* brPhiInit = matlabToBrArray_real(mxGetPr(prhs[0]), dimPhiInit, 3);
	
	// rhoInit
	const mwSize* dimRhoInit = mxGetDimensions(prhs[1]);
	Real* brRhoInit = matlabToBrArray_real(mxGetPr(prhs[1]), dimRhoInit, 3);

	// piInit
	Real* piInit = mxGetPr(prhs[2]);

	// gamInit
	const mwSize* dimGamInit = mxGetDimensions(prhs[3]);
	Real* brGamInit = matlabToBrArray_real(mxGetPr(prhs[3]), dimGamInit, 2);

	// min number of iterations
	Int mini = ((Int*)mxGetData(prhs[4]))[0];

	// max number of iterations
	Int maxi = ((Int*)mxGetData(prhs[5]))[0];

	// gamma
	Real gamma = mxGetScalar(prhs[6]);

	// lambda
	Real lambda = mxGetScalar(prhs[7]);

	// X
	const mwSize* dimX = mxGetDimensions(prhs[8]);
	Real* brX = matlabToBrArray_real(mxGetPr(prhs[8]), dimX, 2);
	
	// Y
	const mwSize* dimY = mxGetDimensions(prhs[9]);
	Real* brY = matlabToBrArray_real(mxGetPr(prhs[9]), dimY, 2);
	
	// tau
	Real tau = mxGetScalar(prhs[10]);

	/////////////
	// OUTPUTS //
	/////////////

	// phi
	const mwSize dimPhi[] = {dimPhiInit[0], dimPhiInit[1], dimPhiInit[2]};
	plhs[0] = mxCreateNumericArray(3,dimPhi,mxDOUBLE_CLASS,mxREAL);
	Real* phi = mxGetPr(plhs[0]);

	// rho
	const mwSize dimRho[] = {dimRhoInit[0], dimRhoInit[1], dimRhoInit[2]};
	plhs[1] = mxCreateNumericArray(3,dimRho,mxDOUBLE_CLASS,mxREAL);
	Real* rho = mxGetPr(plhs[1]);

	// pi
	plhs[2] = mxCreateNumericMatrix(k,1,mxDOUBLE_CLASS,mxREAL);
	Real* pi = mxGetPr(plhs[2]);

	// LLF
	plhs[3] = mxCreateNumericMatrix(maxi,1,mxDOUBLE_CLASS,mxREAL);
	Real* LLF = mxGetPr(plhs[3]);

	// S
	const mwSize dimS[] = {p, m, k};
    plhs[4] = mxCreateNumericArray(3,dimS,mxDOUBLE_CLASS,mxREAL);
	Real* S = mxGetPr(plhs[4]);

	////////////////////
	// Call to EMGLLF //
	////////////////////

	EMGLLF(brPhiInit,brRhoInit,piInit,brGamInit,mini,maxi,gamma,lambda,brX,brY,tau,
		phi,rho,pi,LLF,S,
		n,p,m,k);
	
	free(brPhiInit);
	free(brRhoInit);
	free(brGamInit);
	free(brX);
	free(brY);







}
