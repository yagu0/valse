#ifndef select_EMGrank_H
#define select_EMGrank_H

#include "ioutils.h"

void EMGrank(
	// IN parameters
	const Real* Pi,
	const Real* Rho,
	Int mini,    
	Int maxi,
	const Real* X,
	const Real* Y,
	Real tau,
	const Int* rank,
	// OUT parameters
	Real* phi,
	Real* LLF,
	// additional size parameters
	mwSize n,       
	mwSize p,
	mwSize m,
	mwSize k);

#endif
