// Check if array == refArray
void compareArray(const char* ID, const void* array, const void* refArray, int size, int isInteger);

void compareArray_real(const char* ID, const void* array, const void* refArray, int size);

void compareArray_int(const char* ID, const void* array, const void* refArray, int size);

// Read array by columns (as in MATLAB) and return by-rows encoding
void* readArray(const char* fileName, int isInteger);

int* readArray_int(const char* fileName);

float* readArray_real(const char* fileName);

int read_int(const char* fileName);

float read_real(const char* fileName);
