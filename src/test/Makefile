CC = gcc
CFLAGS = -g -std=gnu99 -Wno-implicit-function-declaration
LDFLAGS = -Lobj/ -lvalse_core
LIB = valse_core.so
LIB_SRC = $(wildcard ../sources/*.c)
LIB_OBJ = $(LIB_SRC:.c=.o)

all: $(LIB) test.EMGLLF test.EMGrank test.constructionModelesLassoMLE test.EMGrank test.constructionModelesLassoRank test.selectionTotale

$(LIB): $(LIB_OBJ)
	$(CC) -o $@ $^

test.EMGLLF: test.EMGLLF.o
	$(CC) -o $@ $^ $(LDFLAGS)

test.constructionModelesLassoMLE: test.constructionModelesLassoMLE.o
	$(CC) -o $@ $^ $(LDFLAGS)

test.EMGrank: test.EMGrank.o
	$(CC) -o $@ $^ $(LDFLAGS)

test.constructionModelesLassoRank: test.constructionModelesLassoRank.o
	$(CC) -o $@ $^ $(LDFLAGS)

test.selectionTotale: test.selectionTotale.o
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -f *.o ../sources/*.o

cclean: clean
	rm -f $(LIB) $(TEST)

.PHONY: all clean cclean
