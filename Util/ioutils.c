#include "ioutils.h"
#include <string.h>
#include <stdio.h>

// Check if array == refArray
void compareArray(const char* ID, const void* array, const void* refArray, mwSize size, int isInteger)
{
	Real EPS = 1e-5; //precision
	printf("Checking %s\n",ID);
	Real maxError = 0.0;
	for (mwSize i=0; i<size; i++)
	{
		Real error = isInteger
			? fabs(((Int*)array)[i] - ((Int*)refArray)[i])
			: fabs(((Real*)array)[i] - ((Real*)refArray)[i]);
		if (error >= maxError)
			maxError = error;
	}
	if (maxError >= EPS)
		printf("    Inaccuracy: max(abs(error)) = %g >= %g\n",maxError,EPS);
	else
		printf("    OK\n");
}

// Next function is a generalization of :
//~ Real* brToMatlabArray(Real* brArray, int dimX, int dimY, int dimZ, int dimW)
//~ {
	//~ Real* matlabArray = (Real*)malloc(dimX*dimY*dimZ*dimW*sizeof(Real));
	//~ for (int u=0; u<dimX*dimY*dimZ*dimW; u++)
	//~ {
		//~ int xIndex = u / (dimY*dimZ*dimW);
		//~ int yIndex = (u % (dimY*dimZ*dimW)) / (dimZ*dimW);
		//~ int zIndex = (u % (dimZ*dimW)) / dimW;
		//~ int wIndex = u % dimW;
		//~ matlabArray[xIndex+yIndex*dimX+zIndex*dimX*dimY+wIndex*dimX*dimY*dimZ] = brArray[u];
	//~ }
	//~ return matlabArray;
//~ }

// Auxiliary to convert from ours ("by-rows") encoding to MATLAB
void* brToMatlabArray(const void* brArray, const mwSize* dimensions, int nbDims, int isInteger)
{
	mwSize totalDim = 1;
	for (int i=0; i<nbDims; i++)
		totalDim *= dimensions[i];
	size_t elementSize = isInteger
		? sizeof(Int)
		: sizeof(Real);
	
	void* matlabArray = malloc(totalDim*elementSize);
	for (mwSize u=0; u<totalDim; u++)
	{
		mwSize prodDimLeft = totalDim;
		mwSize prodDimRight = totalDim / dimensions[0];
		mwSize prodDimInIndex = 1;
		mwSize index = 0;
		for (int v=0; v<nbDims; v++)
		{
			index += ((u % prodDimLeft) / prodDimRight) * prodDimInIndex;
			prodDimInIndex *= dimensions[v];
			prodDimLeft /= dimensions[v];
			if (v+1 < nbDims)
				prodDimRight /= dimensions[v+1];
		}
		if (isInteger)
			((Int*)matlabArray)[index] = ((Int*)brArray)[u];
		else
			((Real*)matlabArray)[index] = ((Real*)brArray)[u];
	}
	return matlabArray;
}

// Next function is a generalization of :
//~ Real* matlabToBrArray(Real* matlabArray, int dimX, int dimY, int dimZ, int dimU)
//~ {
	//~ Real* brArray = (Real*)malloc(dimX*dimY*dimZ*dimU*sizeof(Real));
	//~ for (int u=0; u<dimX*dimY*dimZ*dimU; u++)
	//~ {
		//~ int xIndex = u % dimX;
		//~ int yIndex = (u % (dimX*dimY)) / dimX;
		//~ int zIndex = (u % (dimX*dimY*dimZ)) / (dimX*dimY);
		//~ int uIndex = u / (dimX*dimY*dimZ);
		//~ brArray[xIndex*dimY*dimZ*dimU+yIndex*dimZ*dimU+zIndex*dimU+uIndex] = matlabArray[u];
	//~ }
	//~ return brArray;
//~ }

// Auxiliary to convert from MATLAB encoding to ours ("by-rows")
void* matlabToBrArray(const void* matlabArray, const mwSize* dimensions, int nbDims, int isInteger)
{
	mwSize totalDim = 1;
	for (int i=0; i<nbDims; i++)
		totalDim *= dimensions[i];
	size_t elementSize = isInteger
		? sizeof(Int)
		: sizeof(Real);
	
	void* brArray = malloc(totalDim*elementSize);
	for (mwSize u=0; u<totalDim; u++)
	{
		mwSize prodDimLeft = dimensions[0];
		mwSize prodDimRight = 1;
		mwSize prodDimInIndex = totalDim / dimensions[0];
		mwSize index = 0;
		for (int v=0; v<nbDims; v++)
		{
			index += ((u % prodDimLeft) / prodDimRight) * prodDimInIndex;
			if (v+1 < nbDims)
			{
				prodDimInIndex /= dimensions[v+1];
				prodDimLeft *= dimensions[v+1];
			}
			prodDimRight *= dimensions[v];
		}
		if (isInteger)
			((Int*)brArray)[index] = ((Int*)matlabArray)[u];
		else
			((Real*)brArray)[index] = ((Real*)matlabArray)[u];
	}
	return brArray;
}

// Read array by columns (as in MATLAB) and return by-rows encoding
void* readArray(const char* fileName, const mwSize* dimensions, int nbDims, int isInteger)
{
	// need to prepend '../data/' (not really nice code...)
	char* fullFileName = (char*)calloc(8+strlen(fileName)+1,sizeof(char));
	strcat(fullFileName, "../data/");
	strcat(fullFileName, fileName);
	FILE* file = fopen(fullFileName, "r");
	free(fullFileName);
	
	mwSize totalDim = 1;
	for (int i=0; i<nbDims; i++)
		totalDim *= dimensions[i];
	size_t elementSize = isInteger
		? sizeof(Int)
		: sizeof(Real);
	
	// read all values, and convert them to by-rows matrices format
	void* matlabArray = malloc(totalDim*elementSize);
	char curChar = ' ';
	char bufferNum[64];
	for (mwSize u=0; u<totalDim; u++)
	{
		// position to next non-separator character
		while (!feof(file) && (curChar==' ' || curChar=='\n' || curChar=='\t' || curChar==','))
			curChar = fgetc(file);
		// read number (as a string)
		int bufferIndex = 0;
		while (!feof(file) && curChar!=' ' && curChar!='\n' && curChar!='\t' && curChar!=',')
		{
			bufferNum[bufferIndex++] = curChar;
			curChar = fgetc(file);
		}
		bufferNum[bufferIndex] = 0;
		// transform string into Real, and store it at appropriate location
		if (isInteger)
			((Int*)matlabArray)[u] = atoi(bufferNum);
		else
			((Real*)matlabArray)[u] = atof(bufferNum);
	}
	fclose(file);
	
	void* brArray = matlabToBrArray(matlabArray, dimensions, nbDims, isInteger);
	free(matlabArray);
	return brArray;
}
