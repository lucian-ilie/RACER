
MAKE=make
CXX=g++ -O3 -Wall -fomit-frame-pointer -fopenmp

RACER: RACER.cpp
	$(CXX) RACER.cpp -o $@

clean:
	rm -f *.o
	rm -f RACER
