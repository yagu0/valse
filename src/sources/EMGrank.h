#ifndef valse_EMGrank_H
#define valse_EMGrank_H

void EMGrank(
	// IN parameters
	const double* Pi,
	const double* Rho,
	int mini,
	int maxi,
	const double* X,
	const double* Y,
	double tau,
	const int* rank,
	// OUT parameters
	double* phi,
	double* LLF,
	// additional size parameters
	int n,
	int p,
	int m,
	int k);

#endif
