CC = gcc
CFLAGS = -g -std=gnu99 -Wno-implicit-function-declaration
LDFLAGS = -lm -lgsl -lcblas
TEST_LDFLAGS = -L. libvalse_core.so
LIB = libvalse_core.so
LIB_SRC = $(wildcard ../pkg/src/sources/*.c)
LIB_OBJ = $(LIB_SRC:.c=.o)
INCLUDES = -I../pkg/src/sources
TESTS = test.EMGLLF test.EMGrank

all: $(LIB) $(TESTS)

$(LIB): $(LIB_OBJ)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

test.EMGLLF: $(LIB) test.EMGLLF.o test_utils.o
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

test.EMGrank: $(LIB) test.EMGrank.o test_utils.o
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

%.o: %.c
	$(CC) -fPIC -o $@ -c $< $(CFLAGS) $(INCLUDES)

clean:
	rm -f *.o ../pkg/src/sources/*.o

cclean: clean
	rm -f *.so ../pkg/src/*.so $(TESTS)

.PHONY: all clean cclean
