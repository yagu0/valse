#ifndef valse_EMGrank_H
#define valse_EMGrank_H

void EMGrank_core(
	// IN parameters
	const float* Pi,
	const float* Rho,
	int mini,
	int maxi,
	const float* X,
	const float* Y,
	float tau,
	const int* rank,
	// OUT parameters
	float* phi,
	float* LLF,
	// additional size parameters
	int n,
	int p,
	int m,
	int k);

#endif
