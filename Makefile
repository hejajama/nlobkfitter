AMPLITUDEINCLUDE=$(NLOBKFITTER_LIB_PATH)/amplitudelib
AMPLITUDEPATH=$(NLOBKFITTER_LIB_PATH)/amplitudelib/build/lib
CUBAPATH=$(NLOBKFITTER_LIB_PATH)/Cuba-4.2/lib/
MINUITINC=$(NLOBKFITTER_LIB_PATH)/Minuit2/include/
MINUITLIBDIR=$(NLOBKFITTER_LIB_PATH)/Minuit2/lib/
CXXFLAGS = `gsl-config --cflags` -std=c++11 -O3 -fopenmp -pedantic -I$(AMPLITUDEINCLUDE) -I src/ -I$(MINUITINC)
LDFLAGS = `gsl-config --libs` -L$(MINUITLIBDIR) -lMinuit2 -L$(CUBAPATH) -lcuba -lm

include filelist.m

all: fit

fit: $(OBJECTS) src/main.o
	g++ $(OBJECTS) src/main.o amplitudelib_v2/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o fit $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

fitex: $(OBJECTS) src/main_switches.o 
	g++ $(OBJECTS) src/main_switches.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o fitex $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

swarmscan: $(OBJECTS) src/swarmscan.o 
	g++ $(OBJECTS) src/swarmscan.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o swarmscan $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

lightqscan: $(OBJECTS) src/lightq_scan.o
	g++ $(OBJECTS) src/lightq_scan.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o lightqscan $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)


plottool: $(OBJECTS) src/tool_plot.o 
	g++ $(OBJECTS) src/tool_plot.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o plottool $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

tool: $(OBJECTS) src/toolF2.o 
	g++ $(OBJECTS) src/toolF2.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o tool $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

subunsubtest: $(OBJECTS) src/sub_unsub_test.o 
	g++ $(OBJECTS) src/sub_unsub_test.o $(AMPLITUDEPATH)/libamplitude.a ${MINUITLIBDIR}/libMinuit2.a -o subunsub $(CXXFLAGS) $(PROFFLAGS) $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(AMPLITUDELIBO) src/main.o src/toolF2.o src/main_switches.o src/tool_plot.o src/swarmscan.o src/lightq_scan.o src/sub_unsub_test.o
	rm -f fit tool plottool swarmscan lightqscan subunsub
