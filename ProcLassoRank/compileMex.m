%compile C code (for MATLAB or Octave)
if exist('octave_config_info')
	setenv('CFLAGS','-O2 -std=gnu99 -fopenmp')
	mkoctfile --mex -DOctave -I../Util EMGrank.c EMGrank_interface.c ../Util/ioutils.c -o EMGrank -lm -lgsl -lgslcblas -lgomp
    mkoctfile --mex -DOctave -I../Util constructionModelesLassoRank.c constructionModelesLassoRank_interface.c EMGrank.c ../Util/ioutils.c -o constructionModelesLassoRank -lm -lgsl -lgslcblas -lgomp
else
	mex CFLAGS="\$CFLAGS -O2 -std=gnu99 -fopenmp" -I../Util EMGrank.c EMGrank_interface.c ../Util/ioutils.c -output EMGrank -lm -lgsl -lgslcblas -lgomp
    mex CFLAGS="\$CFLAGS -O2 -std=gnu99 -fopenmp" -I../Util constructionModelesLassoRank.c constructionModelesLassoRank_interface.c EMGrank.c ../Util/ioutils.c -output constructionModelesLassoRank -lm -lgsl -lgslcblas -lgomp
end
