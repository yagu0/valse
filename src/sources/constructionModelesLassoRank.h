#ifndef valse_constructionModelesLassoRank_H
#define valse_constructionModelesLassoRank_H

#include "utils.h"

// Main job on raw inputs (after transformation from mxArray)
void constructionModelesLassoRank_core(
	// IN parameters
	const Real* Pi,
	const Real* Rho,
	int mini,
	int maxi,
	const Real* X,
	const Real* Y,
	Real tau,
	const int* A1,
	int rangmin,
	int rangmax,
	// OUT parameters
	Real* phi,
	Real* lvraisemblance,
	// additional size parameters
	int n,
	int p,
	int m,
	int k,
	int L);

#endif
