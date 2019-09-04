AMPLITUDEINCLUDE=../libs/amplitudelib
AMPLITUDEPATH=../libs/amplitudelib/build/lib
CUBAPATH=../libs/Cuba-4.2/lib/
MINUITINC=../libs/Minuit2/include/
MINUITLIBDIR=/home/henri/Work/libs/Minuit2/lib/
#GSLLIBDIR=/usr/local/gsl-2.3/lib/
CXXFLAGS = `gsl-config --cflags` -std=c++11 -O3 -fopenmp -pedantic -I$(AMPLITUDEINCLUDE) -I src/ -I$(MINUITINC)
LDFLAGS = `gsl-config --libs` -L$(MINUITLIBDIR) -lMinuit2 -L$(CUBAPATH) -lcuba -lm
#PROFFLAGS = -fno-omit-frame-pointer
LIBS_PATH=/home/henri/Work/libs/Minuit2/lib/

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

subunsubtest: $(OBJECTS) src/sub_unsub_test.o
	LD_LIBRARY_PATH=$(LIBS_PATH) g++ $(OBJECTS) src/sub_unsub_test.o $(AMPLITUDEPATH)/libamplitude.a -o subunsub $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(AMPLITUDELIBO) src/main.o src/toolF2.o src/main_switches.o src/tool_plot.o src/sub_unsub_test.o
	rm -f fit tool plottool swarmscan subunsub
