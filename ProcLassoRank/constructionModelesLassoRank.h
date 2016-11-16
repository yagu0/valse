#ifndef select_constructionModelesLassoRank_H
#define select_constructionModelesLassoRank_H

#include "ioutils.h"

// Main job on raw inputs (after transformation from mxArray)
void constructionModelesLassoRank(
	// IN parameters 
	const Real* Pi,
	const Real* Rho,
	Int mini,
	Int maxi,
	const Real* X,
	const Real* Y,
	Real tau,
	const Int* A1,
	Int rangmin,
	Int rangmax,
	// OUT parameters
	Real* phi,
    Real* lvraisemblance,
    // additional size parameters
	mwSize n,         
	mwSize p,
	mwSize m,
	mwSize k,
	mwSize L);

#endif
