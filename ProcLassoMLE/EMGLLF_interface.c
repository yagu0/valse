#include "ioutils.h"
#include "EMGLLF.h"
#include <mex.h>

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
	if (nrhs!=11) 
		mexErrMsgIdAndTxt("select:EMGLLF:nrhs","11 inputs required.");
	if (nlhs!=5) 
		mexErrMsgIdAndTxt("select:EMGLLF:nlhs","5 outputs required.");
	
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
	
	//post-processing: convert by-rows outputs to MATLAB matrices
	Real* mlPhi = brToMatlabArray_real(phi, dimPhi, 3);
	copyArray(mlPhi, phi, dimPhi[0]*dimPhi[1]*dimPhi[2]);
	free(mlPhi);
	
	Real* mlRho = brToMatlabArray_real(rho, dimRho, 3);
	copyArray(mlRho, rho, dimRho[0]*dimRho[1]*dimRho[2]);
	free(mlRho);
    
    Real* mlS = brToMatlabArray_real(S, dimS, 3);
	copyArray(mlS, S, dimS[0]*dimS[1]*dimS[2]);
	free(mlS);
}
