CXX=g++
CFLAGS =-I.
ROOTCFLAGS := $(shell  $(ROOTSYS)/bin/root-config --cflags --libs)
OBJ = AverageAngCorResults.cpp


AverageAngCorResults :  $(OBJ)
	$(CXX) -o $@ $^ $(ROOTCFLAGS) $(CFLAGS)

clean:
	rm -f *.o $(objects) AverageAngCorResults