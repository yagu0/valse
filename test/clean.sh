#!/bin/sh

# Erase roxygen2 generated files
find ../pkg/man/ -type f ! -name 'valse-package.Rd' -delete
rm -f ../pkg/NAMESPACE

# Erase object and library files
rm -f ../pkg/src/*.so
rm -f ../pkg/src/adapters/*.o

make cclean
