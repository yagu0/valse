#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "utils.h"

// Check if array == refArray
void compareArray(const char* ID, const void* array, const void* refArray, int size,
	int isinteger)
{
	Real EPS = 1e-5; //precision
	printf("Checking %s\n",ID);
	Real maxError = 0.0;
	for (int i=0; i<size; i++)
	{
		Real error = isinteger
			? fabs(((int*)array)[i] - ((int*)refArray)[i])
			: fabs(((Real*)array)[i] - ((Real*)refArray)[i]);
		if (error >= maxError)
			maxError = error;
	}
	if (maxError >= EPS)
		printf("    Inaccuracy: max(abs(error)) = %g >= %g\n",maxError,EPS);
	else
		printf("    OK\n");
}

void compareArray_real(const char* ID, const void* array, const void* refArray, int size)
{
	return compareArray(ID, array, refArray, size, 0);
}

void compareArray_int(const char* ID, const void* array, const void* refArray, int size)
{
	return compareArray(ID, array, refArray, size, 1);
}

// Read array by columns (as in MATLAB) and return by-rows encoding
void* readArray(const char* fileName, int isinteger)
{
	// need to prepend 'data/' (not really nice code...)
	char* fullFileName = (char*)calloc(5+strlen(fileName)+1, sizeof(char));
	strcat(fullFileName, "data/");
	strcat(fullFileName, fileName);

	// first pass to know how many elements to allocate
	char* command = (char*)calloc(12+strlen(fullFileName)+8+1, sizeof(char));
	strcat(command, "wc -l ");
	strcat(command, fullFileName);
	FILE *arraySize = popen(command, "r");
	free(command);
	char* bufferNum = (char*)calloc(64, sizeof(char));
	fgets(bufferNum, sizeof(bufferNum), arraySize);
	int n = atoi(bufferNum);
	pclose(arraySize);

	// open file for reading
	FILE* arrayFile = fopen(fullFileName, "r");
	free(fullFileName);

	// read all values, and convert them to by-rows matrices format
	size_t elementSize = isinteger ? sizeof(int) : sizeof(Real);
	void* array = malloc(n*elementSize);
	for (int i=0; i<n; i++)
	{
		fgets(bufferNum, 64, arrayFile);
		// transform buffer content into Real or int, and store it at appropriate location
		if (isinteger)
			((int*)array)[i] = atoi(bufferNum);
		else
			((Real*)array)[i] = atof(bufferNum);
	}
	fclose(arrayFile);
	free(bufferNum);

	return array;
}

int* readArray_int(const char* fileName)
{
	return (int*)readArray(fileName, 1);
}

Real* readArray_real(const char* fileName)
{
	return (Real*)readArray(fileName, 0);
}

int read_int(const char* fileName)
{
	int* array = readArray_int(fileName);
	int res = array[0];
	free(array);
	return res;
}

Real read_real(const char* fileName)
{
	Real* array = readArray_real(fileName);
	Real res = array[0];
	free(array);
	return res;
}
