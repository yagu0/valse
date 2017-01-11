// Check if array == refArray
void compareArray(const char* ID, const void* array, const void* refArray, int size,
	int isInteger)
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

void compareArray_real(const char* ID, const void* array, const void* refArray, int size)
{
	return compareArray(ID, array, refArray, size, 0);
}

void compareArray_int(const char* ID, const void* array, const void* refArray, int size)
{
	return compareArray(ID, array, refArray, size, 1);
}

// Read array by columns (as in MATLAB) and return by-rows encoding
void* readArray(const char* fileName, int isInteger)
{
	// need to prepend '../data/' (not really nice code...)
	char* fullFileName = (char*)calloc(5+strlen(fileName)+1,sizeof(char));
	strcat(fullFileName, "data/");
	strcat(fullFileName, fileName);
	FILE* file = fopen(fullFileName, "r");
	free(fullFileName);

	// first pass to know how many elements to allocate
	// /////................... TODO

	int d = 1;
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

