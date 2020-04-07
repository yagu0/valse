#!/bin/sh
set -e

#Testing procedure for EMGLLF (inside this folder):

algo=$1 #EMGLLF or EMGrank,
        #second arg indicate if rebuild or rebuild+clean requested

if [ "$2" == 'c' ]; then
	#0.1) Clean package + C testing code
	find ../pkg/man/ -type f ! -name 'valse-package.Rd' -delete
	rm -f ../pkg/NAMESPACE
	# Erase object and library files
	rm -f ../pkg/src/*.so
	rm -f ../pkg/src/adapters/*.o
	make cclean
fi

if [ "$2" == 'r' ] || [ "$2" == 'c' ]; then
	#0.2) Install current version of the package (WARNING: roxygen 2 v5.0.1)
	#   --> devtools::install_github('klutometis/roxygen@v5.0.1')
	echo "setwd('../pkg');library(roxygen2);roxygenize('.')" | R --slave
	R CMD INSTALL ../pkg
fi

#1) Generate data using R versions of EMGLLF/EMGrank (slow, but trusted)
echo -e "source('generateRunSaveTest_$algo.R');\n \
		# I'm happy with default values - feel free to give args\n \
		generateRunSaveTest_$algo() " \
	| R --slave

#2) Compile test C code into an executable named "test.$algo"
make test.$algo

#3) Run it with valgrind!
#valgrind
./test.$algo
