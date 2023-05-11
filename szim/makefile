all: szim geom

szim: szim.cpp
	mpicxx  -O3 -std=c++11 -I/opt/mfem/include -I/opt/gslib/include -fopenmp szim.cpp -o szim -L/opt/mfem/lib -lmfem -lrt -llapack -lblas -L/opt/gslib/lib -lgs -lcblas

geom: geom.cpp
	mpicxx  -O3 -std=c++11 -I/opt/mfem/include -fopenmp geom.cpp -o geom -L/opt/mfem/lib -lmfem -lrt -llapack -lblas

