#ifndef select_ioutils_H
#define select_ioutils_H

// Fill an array with zeros
#define zeroArray(array, size)\
{\
	for (int u=0; u<size; u++)\
		array[u] = 0;\
}

// Copy an 1D array
#define copyArray(array, copy, size)\
{\
	for (int u=0; u<size; u++)\
		copy[u] = array[u];\
}

// Matrix Index ; TODO? ncol unused
#define mi(i,j,nrow,ncol)\
	j*nrow + i

// Array Index ; TODO? d3 unused
#define ai(i,j,k,d1,d2,d3)\
	k*d1*d2 + j*d1 + i

#endif
