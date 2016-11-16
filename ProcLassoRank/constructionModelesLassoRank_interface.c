#include "ioutils.h"
#include "constructionModelesLassoRank.h"
#include <mex.h>

#include <stdio.h>

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
	if (nrhs!=10) 
		mexErrMsgIdAndTxt("select:constructionModelesLassoRank:nrhs","10 inputs required.");
	if (nlhs!=2) 
		mexErrMsgIdAndTxt("select:constructionModelesLassoRank:nlhs","3 outputs required.");

	// Get matrices dimensions, to be given to main routine above
	const mwSize n = mxGetDimensions(prhs[4])[0];
	const mwSize p = mxGetDimensions(prhs[4])[1];
	const mwSize m = mxGetDimensions(prhs[1])[0];
	const mwSize k = mxGetDimensions(prhs[1])[2];
	const mwSize L = mxGetDimensions(prhs[7])[1];

	////////////
	// INPUTS //
	////////////

	// pi
	const mwSize* dimPi = mxGetDimensions(prhs[0]);
	Real* brPi = matlabToBrArray_real(mxGetPr(prhs[0]), dimPi, 2);

	// rho
	const mwSize* dimRho = mxGetDimensions(prhs[1]);
	Real* brRho = matlabToBrArray_real(mxGetPr(prhs[1]), dimRho, 4);

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
	
	// A1
	const mwSize* dimA = mxGetDimensions(prhs[7]);
	Int* brA1 = matlabToBrArray_int(mxGetData(prhs[7]), dimA, 2);
	
    //rangmin
    Int rangmin = ((Int*)mxGetData(prhs[8]))[0];
    
    //rangmax
    Int rangmax = ((Int*)mxGetData(prhs[9]))[0];
	
	/////////////
	// OUTPUTS //
	/////////////
	
	// phi
	mwSize Size = pow(rangmax-rangmin+1,k);
	const mwSize dimPhi[] = {p, m, k, L*Size};
	plhs[0] = mxCreateNumericArray(4,dimPhi,mxDOUBLE_CLASS,mxREAL);
	Real* phi = mxGetPr(plhs[0]);

	// lvraisemblance
	const mwSize dimLvraisemblance[] = {L*Size, 2};
	plhs[1] = mxCreateNumericMatrix(dimLvraisemblance[0],dimLvraisemblance[1],mxDOUBLE_CLASS,mxREAL);
	Real* lvraisemblance = mxGetPr(plhs[1]);

	
	//////////////////////////////////////////
	// Call to constructionModelesLassoRank //
	//////////////////////////////////////////

	constructionModelesLassoRank(
        brPi,brRho,mini,maxi,brX,brY,tau,brA1,rangmin,rangmax,
        phi,lvraisemblance,
        n,p,m,k,L);
	
	free(brPi);
	free(brRho);
	free(brX);
	free(brY);
	free(brA1);

	//post-processing: convert by-rows outputs to MATLAB matrices
	Real* mlPhi = brToMatlabArray_real(phi, dimPhi, 4);
	copyArray(mlPhi, phi, dimPhi[0]*dimPhi[1]*dimPhi[2]*dimPhi[3]);
	free(mlPhi);
	
	Real* mlLvraisemblance = brToMatlabArray_real(lvraisemblance, dimLvraisemblance, 2);
	copyArray(mlLvraisemblance, lvraisemblance, dimLvraisemblance[0]*dimLvraisemblance[1]);
	free(mlLvraisemblance);

}
