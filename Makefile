CXX=g++
CXXFLAGS=-O3 -fopenmp
DEF= -D___UNIFORM


0d: 0d.o ap.o cell.o recsubcell.o subcell.o log.o
	$(CXX) -o $@ $(CXXFLAGS) 0d.o ap.o cell.o recsubcell.o subcell.o log.o

.cc.o:
	$(CXX) -c -MMD -MP $< $(CXXFLAGS)

clean:
	rm *.o
