#ifndef valse_constructionModelesLassoRank_H
#define valse_constructionModelesLassoRank_H

// Main job on raw inputs (after transformation from mxArray)
void constructionModelesLassoRank_core(
	// IN parameters
	const double* Pi,
	const double* Rho,
	int mini,
	int maxi,
	const double* X,
	const double* Y,
	double tau,
	const int* A1,
	int rangmin,
	int rangmax,
	// OUT parameters
	double* phi,
	double* lvraisemblance,
	// additional size parameters
	int n,
	int p,
	int m,
	int k,
	int L);

#endif
