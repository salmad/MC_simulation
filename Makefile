CC=mpic++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp src/molecule.cpp src/polymer.cpp src/sim_system.cpp src/star.cpp integrate.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mc_star

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@