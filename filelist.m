SOURCES = src/dipole.cpp src/data.cpp src/solver.cpp src/ic.cpp src/mv.cpp src/ic_datafile.cpp src/nlobk_config.cpp src/nlodis_config.cpp src/nlodissigmar.cpp src/nlodissigmar_massiveq.cpp src/helper.cpp
OBJECTS=$(SOURCES:.cpp=.o)
DEPENDENCIES = $(SOURCES:.cpp=.hpp)
