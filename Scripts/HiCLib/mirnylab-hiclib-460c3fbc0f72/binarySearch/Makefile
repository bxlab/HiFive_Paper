all:
	cython --cplus fastBinSearch.pyx  
	g++  -g -Ofast  -fpic -c -std=c++0x -O2    mycode.cpp -o mycode.o -fopenmp
	g++  -g -Ofast  -fpic -c -std=c++0x  fastBinSearch.cpp -o fastBinSearch.o `python-config --includes` 
	g++ -g -Ofast  -std=c++0x  -shared -o fastBinSearch.so fastBinSearch.o mycode.o `python-config --libs` -lgomp
	cp fastBinSearch.so ../src/hiclib

clean:
	rm -f mycodepy.cpp *.o *.so
