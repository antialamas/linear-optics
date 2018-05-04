CC = g++
CFLAGS = -Ofast -funroll-loops -c
LFLAGS = -Ofast -funroll-loops
OBJS = LinearOpticalTransform.o main.o GrayCode.o

all: LinearOpticalSimulation

LinearOpticalSimulation: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o LinearOpticalSimulation

GrayCode.o: GrayCode.cpp
	$(CC) $(CFLAGS) GrayCode.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

LinearOpticalTransform.o: LinearOpticalTransform.cpp
	$(CC) $(CFLAGS) LinearOpticalTransform.cpp

clean:
	rm *.o LinearOpticalSimulation
