#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// Check if array == refArray
void compareArray(const char* ID, const void* array, const void* refArray, int size,
	int isinteger)
{
	float EPS = 1e-5; //precision
	printf("Checking %s\n",ID);
	float maxError = 0.0;
	for (int i=0; i<size; i++)
	{
		float error = isinteger
			? fabs(((int*)array)[i] - ((int*)refArray)[i])
			: fabs(((float*)array)[i] - ((float*)refArray)[i]);
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
	// need to prepend '../data/' (not really nice code...)
	char* fullFileName = (char*)calloc(5+strlen(fileName)+1, sizeof(char));
	strcat(fullFileName, "data/");
	strcat(fullFileName, fileName);

	// first pass to know how many elements to allocate
	char* command = (char*)calloc(12+strlen(fullFileName)+8+1, sizeof(char));
	strcat(command, "grep -o ' ' ");
	strcat(command, fullFileName);
	strcat(command, " | wc -l");
	FILE *countSpaces = popen(command, "r");
	char* buffer = (char*)calloc(32, sizeof(char));
	fgets(buffer, sizeof(buffer), countSpaces);
	int n = atoi(buffer) + 1;
	free(buffer);
	pclose(countSpaces);

	// open file for reading
	FILE* file = fopen(fullFileName, "r");
	free(fullFileName);

	int d = 1;
	size_t elementSize = isinteger
		? sizeof(int)
		: sizeof(float);

	// read all values, and convert them to by-rows matrices format
	void* array = malloc(n*elementSize);
	char curChar = ' ';
	char bufferNum[64];
	for (int u=0; u<n; u++)
	{
		// read number (as a string)
		int bufferIndex = 0;
		while (!feof(file) && curChar!=' ')
		{
			curChar = fgetc(file);
			bufferNum[bufferIndex++] = curChar;
		}
		bufferNum[bufferIndex] = 0;
		// transform string into float, and store it at appropriate location
		if (isinteger)
			((int*)array)[u] = atoi(bufferNum);
		else
			((float*)array)[u] = atof(bufferNum);
		// position to next non-separator character
		curChar = fgetc(file);
	}
	fclose(file);

	return array;
}

int* readArray_int(const char* fileName)
{
	return (int*)readArray(fileName, 1);
}

float* readArray_real(const char* fileName)
{
	return (float*)readArray(fileName, 0);
}

int read_int(const char* fileName)
{
	return readArray_int(fileName)[0];
}

float read_real(const char* fileName)
{
	return readArray_real(fileName)[0];
}
