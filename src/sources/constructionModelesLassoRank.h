#ifndef valse_constructionModelesLassoRank_H
#define valse_constructionModelesLassoRank_H

// Main job on raw inputs (after transformation from mxArray)
void constructionModelesLassoRank_core(
	// IN parameters
	const float* Pi,
	const float* Rho,
	int mini,
	int maxi,
	const float* X,
	const float* Y,
	float tau,
	const int* A1,
	int rangmin,
	int rangmax,
	// OUT parameters
	float* phi,
	float* lvraisemblance,
	// additional size parameters
	int n,
	int p,
	int m,
	int k,
	int L);

#endif
