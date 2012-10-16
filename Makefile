CXX=g++
CXXFLAGS = -pg -I. -I/opt/local/include/ -I/usr/local/include/quakelib-1.0
LDFLAGS = -L/opt/local/lib/ -lgsl -lsundials_cvode -lsundials_nvecserial
DEPS = BlockData.h Equations.h RateState.h Solver.h Spline.h
OBJ = Equations.o RateState.o Solver.o Spline.o main.o

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

vicars: $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -rf *.o vicars


