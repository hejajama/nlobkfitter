AMPLITUDEINCLUDE=/nashome1/hejohann/Libs/amplitudelib
AMPLITUDEPATH=/nashome1/hejohann/Libs/amplitudelib/build/lib
CUBAPATH=/nashome1/hejohann/Libs/Cuba-4.2/
MINUITINC=/nashome1/hejohann/Libs/Minuit2/include/
MINUITLIBDIR=/nashome1/hejohann/Libs/Minuit2/lib
GSLLIBDIR=#/usr/local/gsl-2.3/lib/
CXXFLAGS = `gsl-config --cflags` -std=c++11 -O3 -fopenmp -I$(AMPLITUDEINCLUDE) -I src/ -I$(MINUITINC)
LDFLAGS = `gsl-config --libs` -L$(MINUITLIBDIR) -lMinuit2 -L$(CUBAPATH) -lcuba -lm
#PROFFLAGS = -fno-omit-frame-pointer
LIBS_PATH=/nashome1/hejohann/Libs/Minuit2/lib

include filelist.m

all: fit

fit: $(OBJECTS) src/main_switches.o 
	LD_LIBRARY_PATH=$(LIBS_PATH) g++ $(OBJECTS) src/main_switches.o $(AMPLITUDEPATH)/libamplitude.a -o fit $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

swarmscan: $(OBJECTS) src/swarmscan.o 
	LD_LIBRARY_PATH=$(LIBS_PATH) g++ $(OBJECTS) src/swarmscan.o $(AMPLITUDEPATH)/libamplitude.a -o swarmscan $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)


plottool: $(OBJECTS) src/tool_plot.o 
	LD_LIBRARY_PATH=$(LIBS_PATH) g++ $(OBJECTS) src/tool_plot.o $(AMPLITUDEPATH)/libamplitude.a -o plottool $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

tool: $(OBJECTS) src/toolF2.o 
	LD_LIBRARY_PATH=$(LIBS_PATH) g++ $(OBJECTS) src/toolF2.o $(AMPLITUDEPATH)/libamplitude.a -o tool $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(AMPLITUDELIBO) src/main.o src/toolF2.o src/main_switches.o src/tool_plot.o
	rm -f fit tool plottool swarmscan
