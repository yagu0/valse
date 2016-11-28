#include "ioutils.h"
#include "EMGrank.h"
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
	if (nrhs!=8) 
		mexErrMsgIdAndTxt("select:EMGrank:nrhs","8 inputs required.");
	if (nlhs!=2) 
		mexErrMsgIdAndTxt("select:EMGrank:nlhs","3 outputs required.");

	// Get matrices dimensions
	const mwSize n = mxGetDimensions(prhs[4])[0];
	const mwSize p = mxGetDimensions(prhs[4])[1];
	const mwSize m = mxGetDimensions(prhs[1])[0];
	const mwSize k = mxGetDimensions(prhs[1])[2];

	////////////
	// INPUTS //
	////////////

	// Pi
	Real* Pi = mxGetPr(prhs[0]);

	// Rho
	const mwSize* dimRho = mxGetDimensions(prhs[1]);
	Real* brRho = matlabToBrArray_real(mxGetPr(prhs[1]), dimRho, 3);
	
	// min number of iterations
	Int mini = ((Int*)mxGetData(prhs[2]))[0];

	// max number of iterations
	Int maxi = ((Int*)mxGetData(prhs[3]))[0];

	// X
	const mwSize* dimX = mxGetDimensions(prhs[4]);
	Real* brX = matlabToBrArray_real(mxGetPr(prhs[4]), dimX, 2);
	
	// Y
	const mwSize* dimY = mxGetDimensions(prhs[5]);
	Real* brY = matlabToBrArray_real(mxGetPr(prhs[5]), dimY, 2);

	// tau
	Real tau = mxGetScalar(prhs[6]);

	// rank
	Int* rank = (Int*)mxGetData(prhs[7]);
	
	/////////////
	// OUTPUTS //
	/////////////

	// phi
	const mwSize dimPhi[] = {p, m, k};
	plhs[0] = mxCreateNumericArray(3,dimPhi,mxDOUBLE_CLASS,mxREAL);
	Real* phi = mxGetPr(plhs[0]);

	// LLF
	plhs[1] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
	Real* LLF = mxGetPr(plhs[1]);

	/////////////////////
	// Call to EMGrank //
	/////////////////////

	EMGrank(Pi,brRho,mini,maxi,brX,brY,tau,rank,
		phi,LLF,
		n,p,m,k);

	free(brRho);
	free(brX);
	free(brY);
	
	//post-processing: convert by-rows outputs to MATLAB matrices
	Real* mlPhi = brToMatlabArray_real(phi, dimPhi, 3);
	copyArray(mlPhi, phi, dimPhi[0]*dimPhi[1]*dimPhi[2]);
	free(mlPhi);
}
