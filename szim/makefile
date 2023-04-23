all: szim geom

szim: szim.cpp
	mpicxx  -O3 -std=c++11 -I/opt/hypre/include -I/opt/metis/include -I/opt/mfem/include -fopenmp szim.cpp -o szim -L/opt/mfem/lib -lmfem -L/opt/hypre/lib -lHYPRE -L/opt/metis/lib -lmetis -lrt

geom: geom.cpp
	mpicxx  -O3 -std=c++11 -I/opt/hypre/include -I/opt/metis/include -I/opt/mfem/include -fopenmp geom.cpp -o geom -L/opt/mfem/lib -lmfem -L/opt/hypre/lib -lHYPRE -L/opt/metis/lib -lmetis -lrt
	
