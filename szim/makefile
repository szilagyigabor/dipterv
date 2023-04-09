all: szim geom

szim: szim.cpp
	mpicxx  -O3 -std=c++11 -I/usr/include -I/usr/local/include -fopenmp szim.cpp -o szim -lmfem -L/usr/lib -lHYPRE -L/usr/local/lib -lmetis -lrt

geom: geom.cpp
	mpicxx  -O3 -std=c++11 -I/usr/include -I/usr/local/include -fopenmp geom.cpp -o geom -lmfem -L/usr/lib -lHYPRE -L/usr/local/lib -lmetis -lrt
	
