#nlobkfitter installation

## Dependencies:
    - cmake
    - GSL
    - amplitudelib_v2
    - Minuit2
    - Cuba (MC integration library)

## Procedure

    - install GSL
    - build amplitudelib_v2 by:
	-- clone https://github.com/hejajama/amplitudelib
	-- mkdir [amplitudelib/]build && cd build && cmake .. && make
    - build Minuit2:
	-- download from CERN homepage
	-- ./configure --disable-openmp --prefix=INSTALLDIR
	-- make -j[N] && make install
    - build Cuba:
	-- download
	-- ./configure --prefix=INSTALLDIR
	-- make && make install
    - build nlobkfitter & tools:
	-- clone https://github.com/hejajama/nlobkfitter
	-- make or cmake
