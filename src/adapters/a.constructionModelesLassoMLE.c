#include "ioutils.h"
#include "constructionModelesLassoMLE.h"
#include <mex.h>

#include <stdio.h>

// nlhs, nrhs: resp. numbers of out and in parameters.
// plhs: array of out parameters, each being a mxArray
// plhs: array of in parameters (immutable), each being a mxArray
//
// MATLAB translates a call [A,B] = fun(C,D) into mexFunction(2,{A,B},2,{C,D}). 
// Then mxArrayS are adapted to be passed to a regular C function, 
// and the results are translated back to mxArrayS into plhs.
void mexFunction(
	int nlhs, 
	mxArray* plhs[],
	int nrhs, 
	const mxArray* prhs[])
{
	// Basic sanity checks
	if (nrhs!=14) 
		mexErrMsgIdAndTxt("select:constructionModelesLassoMLE:nrhs","14 inputs required.");
	if (nlhs!=4) 
		mexErrMsgIdAndTxt("select:constructionModelesLassoMLE:nlhs","4 outputs required.");

	// Get matrices dimensions
	const mwSize n = mxGetDimensions(prhs[8])[0];
	const mwSize p = mxGetDimensions(prhs[0])[0];
	const mwSize m = mxGetDimensions(prhs[0])[1];
	const mwSize k = mxGetDimensions(prhs[0])[2];
	const mwSize L = mxGetNumberOfElements(prhs[7]);
	
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

	// glambda
	Real* glambda = mxGetPr(prhs[7]);

	// X
	const mwSize* dimX = mxGetDimensions(prhs[8]);
	Real* brX = matlabToBrArray_real(mxGetPr(prhs[8]), dimX, 2);

	// Y
	const mwSize* dimY = mxGetDimensions(prhs[9]);
	Real* brY = matlabToBrArray_real(mxGetPr(prhs[9]), dimY, 2);
	
	//seuil
	Real seuil = mxGetScalar(prhs[10]);

	// tau
	Real tau = mxGetScalar(prhs[11]);
	
	// A1
	const mwSize* dimA = mxGetDimensions(prhs[12]);
	Int* brA1 = matlabToBrArray_int(mxGetData(prhs[12]), dimA, 3);
	
	// A2
	Int* brA2 = matlabToBrArray_int(mxGetData(prhs[13]), dimA, 3);
	
	/////////////
	// OUTPUTS //
	/////////////
	
	// phi
	const mwSize dimPhi[] = {dimPhiInit[0], dimPhiInit[1], dimPhiInit[2], L};
	plhs[0] = mxCreateNumericArray(4,dimPhi,mxDOUBLE_CLASS,mxREAL);
	Real* phi = mxGetPr(plhs[0]);
	
	// rho
	const mwSize dimRho[] = {dimRhoInit[0], dimRhoInit[1], dimRhoInit[2], L};
	plhs[1] = mxCreateNumericArray(4,dimRho,mxDOUBLE_CLASS,mxREAL);
	Real* rho = mxGetPr(plhs[1]);

	// pi
	const mwSize dimPi[] = {k, L};
	plhs[2] = mxCreateNumericMatrix(dimPi[0],dimPi[1],mxDOUBLE_CLASS,mxREAL);
	Real* pi = mxGetPr(plhs[2]);

	// lvraisemblance
	const mwSize dimLvraisemblance[] = {L, 2};
	plhs[3] = mxCreateNumericMatrix(L, 2, mxDOUBLE_CLASS,mxREAL);
	Real* lvraisemblance = mxGetPr(plhs[3]);

	/////////////////////////////////////////
	// Call to constructionModelesLassoMLE //
	/////////////////////////////////////////

	constructionModelesLassoMLE(
		brPhiInit,brRhoInit,piInit,brGamInit,mini,maxi,gamma,glambda,brX,brY,seuil,tau,brA1,brA2,
		phi,rho,pi,lvraisemblance,
		n,p,m,k,L);
	
	free(brPhiInit);
	free(brRhoInit);
	free(brGamInit);
	free(brX);
	free(brY);
	free(brA1);
	free(brA2);
	
	//post-processing: convert by-rows outputs to MATLAB matrices
	Real* mlPhi = brToMatlabArray_real(phi, dimPhi, 4);
	copyArray(mlPhi, phi, dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]);
	free(mlPhi);
	
	Real* mlRho = brToMatlabArray_real(rho, dimRho, 4);
	copyArray(mlRho, rho, dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]);
	free(mlRho);

	Real* mlPi = brToMatlabArray_real(pi, dimPi, 2);
	copyArray(mlPi, pi, dimPi[0]*dimPi[1]);
	free(mlPi);

	Real* mlLvraisemblance = brToMatlabArray_real(lvraisemblance, dimLvraisemblance, 2);
	copyArray(mlLvraisemblance, lvraisemblance, dimLvraisemblance[0]*dimLvraisemblance[1]);
	free(mlLvraisemblance);
}
