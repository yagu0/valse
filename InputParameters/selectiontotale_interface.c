#include "ioutils.h"
#include "selectiontotale.h"
#include "EMGLLF.h"
#include <mex.h>

// nlhs, nrhs: resp. numbers of out and in parameters.
// plhs: array of out parameters, each being a mxArray
// plhs: array of in parameters (immutable), each being a mxArray
//
// MATLAB translates a call [A,B] = fun(C,D) into mexFunction(2,{A,B},2,{C,D}). 
// Then mxArrayS are adapted to be passed to a regular C function, 
// and the results are translated back to mxArrayS into plhs.
void mexFunction(int nlhs, mxArray* plhs[], 
                 int nrhs, const mxArray* prhs[])
{
	// Basic sanity checks
	if (nrhs!=12) 
		mexErrMsgIdAndTxt("select:selectiontotale:nrhs","12 inputs required.");
	if (nlhs!=4) 
		mexErrMsgIdAndTxt("select:selectiontotale:nlhs","4 outputs required.");

	// Get matrices dimensions, to be given to main routine above
	const mwSize n = mxGetDimensions(prhs[8])[0];
	const mwSize p = mxGetDimensions(prhs[8])[1];
	const mwSize m = mxGetDimensions(prhs[1])[0];
	const mwSize k = mxGetDimensions(prhs[1])[2];
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

	/////////////
	// OUTPUTS //
	/////////////
	
	// A1
	mwSize dimA[] = {p,m+1,L};
	plhs[0] = mxCreateNumericArray(3,dimA,mxGetClassID(prhs[4]),mxREAL);
	Int* A1 = (Int*)mxGetData(plhs[0]);
	
	// A2
	plhs[1] = mxCreateNumericArray(3,dimA,mxGetClassID(prhs[4]),mxREAL);
	Int* A2 = (Int*)mxGetData(plhs[1]);
	
	// rho
	const mwSize dimRho[] = {dimRhoInit[0], dimRhoInit[1], dimRhoInit[2], L};
	plhs[2] = mxCreateNumericArray(4,dimRho,mxDOUBLE_CLASS,mxREAL);
	Real* Rho = mxGetPr(plhs[2]);

	// pi
	const mwSize dimPi[] = {k, L};
	plhs[3] = mxCreateNumericMatrix(dimPi[0],dimPi[1],mxDOUBLE_CLASS,mxREAL);
	double* Pi = mxGetPr(plhs[3]);

	/////////////////////////////
	// Call to selectiontotale //
	/////////////////////////////

	selectiontotale(brPhiInit,brRhoInit,piInit,brGamInit,mini,maxi,gamma,glambda,brX,brY,seuil,tau,
		A1,A2,Rho,Pi,
		n,p,m,k,L);
	
	free(brPhiInit);
	free(brRhoInit);
	free(brGamInit);
	free(brX);
	free(brY);
	
	//post-processing: convert by-rows outputs to MATLAB matrices
	Int* mlA1 = brToMatlabArray_int(A1,dimA,3);
	copyArray(mlA1,A1,dimA[0]*dimA[1]*dimA[2]);
	free(mlA1);
	
	Int* mlA2 = brToMatlabArray_int(A2,dimA,3);
	copyArray(mlA2,A2,dimA[0]*dimA[1]*dimA[2]);
	free(mlA2);
	
	Real* mlRho = brToMatlabArray_real(Rho, dimRho, 4);
	copyArray(mlRho, Rho, dimRho[0]*dimRho[1]*dimRho[2]*dimRho[3]);
	free(mlRho);

	Real* mlPi = brToMatlabArray_real(Pi, dimPi, 2);
	copyArray(mlPi, Pi, dimPi[0]*dimPi[1]);
	free(mlPi);
}
