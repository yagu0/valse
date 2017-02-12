#!/bin/sh

rm -f src/*.so
rm -f src/adapters/*.o
rm -f src/sources/*.o
cd src/test && make cclean && cd ../..
