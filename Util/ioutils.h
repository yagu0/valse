#ifndef select_ioutils_H
#define select_ioutils_H

#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <uchar.h> //for type wchar16_t

// Include header for mwSize type
#ifdef Octave
#include <mex.h>
#else
#include <tmwtypes.h>
#endif

// CHUNK_SIZE = number of lambda values to be treated sequentially by a single core
#define CHUNK_SIZE 1

// integer type chosen in MATLAB (to be tuned: 32 bits should be enough)
typedef int64_t Int;

// real number type chosen in MATLAB (default: double)
typedef double Real;

#ifndef M_PI
#define M_PI 3.141592653589793 
#endif

// Fill an array with zeros
#define zeroArray(array, size)\
{\
	for (Int u=0; u<size; u++)\
		array[u] = 0;\
}

// Copy an 1D array
#define copyArray(array, copy, size)\
{\
	for (Int u=0; u<size; u++)\
		copy[u] = array[u];\
}

// Check if array == refArray
void compareArray(const char* ID, const void* array, const void* refArray, mwSize size, int isInteger);

#define compareArray_int(ID, array, refArray, size)\
	compareArray(ID, array, refArray, size, 1)
#define compareArray_real(ID, array, refArray, size)\
	compareArray(ID, array, refArray, size, 0)

// Auxiliary to convert from ours ("by-rows") encoding to MATLAB
void* brToMatlabArray(const void* brArray, const mwSize* dimensions, int nbDims, int isInteger);

#define brToMatlabArray_int(brArray, dimensions, nbDims)\
	(Int*)brToMatlabArray(brArray, dimensions, nbDims, 1)
#define brToMatlabArray_real(brArray, dimensions, nbDims)\
	(Real*)brToMatlabArray(brArray, dimensions, nbDims, 0)

// Auxiliary to convert from MATLAB encoding to ours ("by-rows")
void* matlabToBrArray(const void* matlabArray, const mwSize* dimensions, int nbDims, int isInteger);

#define matlabToBrArray_int(matlabArray, dimensions, nbDims)\
	(Int*)matlabToBrArray(matlabArray, dimensions, nbDims, 1)
#define matlabToBrArray_real(matlabArray, dimensions, nbDims)\
	(Real*)matlabToBrArray(matlabArray, dimensions, nbDims, 0)

// Read array by columns (as in MATLAB) and return by-rows encoding
void* readArray(const char* fileName, const mwSize* dimensions, int nbDims, int isInteger);

#define readArray_int(fileName, dimensions, nbDims)\
	(Int*)readArray(fileName, dimensions, nbDims, 1)
#define readArray_real(fileName, dimensions, nbDims)\
	(Real*)readArray(fileName, dimensions, nbDims, 0)

#endif
