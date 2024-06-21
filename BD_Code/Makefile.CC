#CC=pgc++
CC=mpiCC
#CC=g++
#CC=mpic++

CFLAGS=-c 
LDFLAGS= -lm -w
SOURCES=  BD_Rouse_UniExtension_MPI_CC.c
OBJECTS=  BD_Rouse_UniExtension_MPI_CC.o
EXECUTABLE=BD_Rouse_UniExtension_MPI_CC

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(LDFLAGS) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE) 
