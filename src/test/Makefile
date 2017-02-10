CC = gcc
CFLAGS = -g -std=gnu99 -Wno-implicit-function-declaration
LDFLAGS = -lm -lgsl -lcblas -lgomp
TEST_LDFLAGS = -L. libvalse_core.so
LIB = libvalse_core.so
LIB_SRC = $(wildcard ../sources/*.c)
LIB_OBJ = $(LIB_SRC:.c=.o)
INCLUDES = -I../sources

all: $(LIB) test.EMGLLF test.EMGrank test.constructionModelesLassoMLE test.EMGrank test.constructionModelesLassoRank test.selectionTotale

$(LIB): $(LIB_OBJ)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

test.EMGLLF: test.EMGLLF.o test_utils.o
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

test.constructionModelesLassoMLE: test.constructionModelesLassoMLE.o test_utils.o
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

test.EMGrank: test.EMGrank.o test_utils.o
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

test.constructionModelesLassoRank: test.constructionModelesLassoRank.o test_utils.o
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

test.selectionTotale: test.selectionTotale.o test_utils.o
	$(CC) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

%.o: %.c
	$(CC) -fPIC -o $@ -c $< $(CFLAGS) $(INCLUDES)

clean:
	rm -f *.o ../sources/*.o

cclean: clean
	rm -f $(LIB) $(TEST)

.PHONY: all clean cclean
