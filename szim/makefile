all: szim geom

#szim: szim.cpp
#	mpicxx  -O3 -std=c++11 -I/usr/include -I/usr/include/hypre -I/usr/local/include -fopenmp szim.cpp -o szim -lmfem -L/usr/lib -lHYPRE -L/usr/local/lib -lmetis -lrt
#
#geom: geom.cpp
#	mpicxx  -O3 -std=c++11 -I/lib64 -I/usr/include -I/usr/include/hypre -I/usr/local/include -L/usr/local/lib -fopenmp geom.cpp -o geom -lmfem -L/usr/lib -lHYPRE -lmetis -lrt
	
szim: szim.cpp
	mpicxx  -O3 -std=c++11 -I/usr/include -I/usr/include/hypre -I/usr/local/include -fopenmp szim.cpp -o szim -lmfem -L/usr/lib -lHYPRE -L/usr/local/lib -lmetis -lrt

geom: geom.cpp
	mpicxx  -O3 -std=c++11 -I/opt/hypre/include -I/opt/metis/include -I/opt/mfem/include -fopenmp geom.cpp -o geom -L/opt/mfem/lib -lmfem -L/opt/hypre/lib -lHYPRE -L/opt/metis/lib -lmetis -lrt
	
