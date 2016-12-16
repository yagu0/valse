#ifndef valse_utils_H
#define valse_utils_H

/*******************
 * tune parallelism
 *******************/

// Number of OpenMP threads
#define OMP_NUM_THREADS 8

// CHUNK_SIZE = number of lambda values to be treated sequentially by a single core
#define CHUNK_SIZE 1

/*******************************
 * Matrix and arrays indexation
 *******************************/

// Matrix Index ; TODO? ncol unused
#define mi(i,j,nrow,ncol)\
	j*nrow + i

// Array Index ; TODO? d3 unused
#define ai(i,j,k,d1,d2,d3)\
	k*d1*d2 + j*d1 + i

// Array4 Index ; TODO? ...
#define ai4(i,j,k,m,d1,d2,d3,d4)\
	m*d1*d2*d3 + k*d1*d2 + j*d1 + i

// Array5 Index ; TODO? ...
#define ai5(i,j,k,m,n,d1,d2,d3,d4,d5)\
	n*d1*d2*d3*d4 + m*d1*d2*d3 + k*d1*d2 + j*d1 + i

/*************************
 * Array copy & "zeroing"
 ************************/

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

#endif
