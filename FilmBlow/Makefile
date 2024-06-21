
#CC=mpiCC
CC=gcc
CFLAGS=-c 
LDFLAGS=-lm  -w  -lstdc++
EXE_name=FILMBLOW
SOURCES= ./RDP_history.cpp ./Rate_RDP_Cost.cpp ./PSO.cpp
OBJECTS=  RDP_history.o  Rate_RDP_Cost.o PSO.o 
EXECUTABLE=$(EXE_name)

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(LDFLAGS) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE)
