SOURCES = src/dipole.cpp src/data.cpp src/solver.cpp src/ic.cpp src/mv.cpp src/mv-nucl.cpp src/ic_datafile.cpp src/nlobk_config.cpp src/nlodis_config.cpp src/helper.cpp
OBJECTS=$(SOURCES:.cpp=.o)
DEPENDENCIES = $(SOURCES:.cpp=.hpp)
