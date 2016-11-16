%compile C code (for MATLAB or Octave)
if exist('octave_config_info')
	setenv('CFLAGS','-O2 -std=gnu99 -fopenmp')
	mkoctfile --mex -DOctave -I../Util EMGLLF.c EMGLLF_interface.c ../Util/ioutils.c -o EMGLLF -lm -lgsl -lgslcblas -lgomp
    mkoctfile --mex -DOctave -I../Util constructionModelesLassoMLE.c constructionModelesLassoMLE_interface.c EMGLLF.c ../Util/ioutils.c -o constructionModelesLassoMLE -lm -lgsl -lgslcblas -lgomp
else
	mex CFLAGS="\$CFLAGS -O2 -std=gnu99 -fopenmp" -I../Util EMGLLF.c EMGLLF_interface.c ../Util/ioutils.c -output EMGLLF -lm -lgsl -lgslcblas -lgomp
    mex CFLAGS="\$CFLAGS -O2 -std=gnu99 -fopenmp" -I../Util constructionModelesLassoMLE.c constructionModelesLassoMLE_interface.c EMGLLF.c ../Util/ioutils.c -output constructionModelesLassoMLE -lm -lgsl -lgslcblas -lgomp
end
