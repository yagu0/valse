%compile C code (for MATLAB or Octave)
if exist('octave_config_info')
	setenv('CFLAGS','-O2 -std=gnu99 -fopenmp')
	mkoctfile --mex -DOctave -I../Util -I../ProcLassoMLE selectiontotale.c selectiontotale_interface.c ../Util/ioutils.c ../ProcLassoMLE/EMGLLF.c -o selectiontotale -lm -lgsl -lgslcblas -lgomp
else
	mex CFLAGS="\$CFLAGS -O2 -std=gnu99 -fopenmp" -I../Util -I../ProcLassoMLE selectiontotale.c selectiontotale_interface.c ../Util/ioutils.c ../ProcLassoMLE/EMGLLF.c -output selectiontotale -lm -lgsl -lgslcblas -lgomp
end
