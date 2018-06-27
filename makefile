
MAKE=make
CXX=g++ -O3 -Wall -fomit-frame-pointer -fopenmp

racer: RACER.cpp
	$(CXX) RACER.cpp -o $@

clean:
	rm -f *.o
	rm -f racer
