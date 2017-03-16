#!/bin/sh

rm -f ../pkg/src/*.so
rm -f ../pkg/src/adapters/*.o
rm -f ../pkg/src/sources/*.o
make cclean
