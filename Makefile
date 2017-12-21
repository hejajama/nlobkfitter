CUBAPATH=../../Cuba-4.2
MINUITINC=../../minuit/include/
MINUITLIBDIR=../../minuit/lib/
CXXFLAGS = `gsl-config --cflags` -std=c++11 -O2 -fopenmp -pedantic -I amplitudelib_v2/ -I src/ -I$(MINUITINC)
LDFLAGS = `gsl-config --libs` -L$(MINUITLIBDIR)  -L$(CUBAPATH) -lcuba -lm
#PROFFLAGS = -fno-omit-frame-pointer

include filelist.m

all: fit

fit: $(OBJECTS) src/main.o 
	g++ $(OBJECTS) src/main.o amplitudelib_v2/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o fit $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(AMPLITUDELIBO) src/main.o
	rm -f fit	
