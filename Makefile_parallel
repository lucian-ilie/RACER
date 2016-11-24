
MAKE=make
CXX=g++ -O3 -Wall -fomit-frame-pointer -fopenmp

RACER: RACER_parallel.cpp
	$(CXX) RACER_parallel.cpp -o $@

clean:
	rm -f *.o
	rm -f RACER
