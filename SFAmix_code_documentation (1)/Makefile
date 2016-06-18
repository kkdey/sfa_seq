CC=g++
CFLAGS=-O2 -march=x86-64 -ffast-math -fopenmp
#CFLAGS=-O2 -march=x86-64 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -g -pg -fopenmp
LINKER=-L/your linker files for gsl/
INCLUDE=-I/folder of the header files/
#LIBS=-lstdc++ -lgfortran -lgsl -lgslcblas -lpthread -lm -static
LIBS=-lgsl -lgslcblas -lm

OBJECTS=SFAmix.o

main: $(OBJECTS)
	g++ $(CFLAGS) $(OBJECTS) -o SFAmix $(LINKER) $(LIBS)
SFAmix.o: SFAmix.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c SFAmix.cpp
clean:
	rm *.o