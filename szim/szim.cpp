#include "mfem.hpp"
#include <stdio.h>

using namespace std;
using namespace mfem;
//using namespace mfem::common;

int main(int argc, char *argv[])
{
    
    Mpi::Init();
    Hypre::Init();
 
    const char *mesh_file = "discontinuity.mesh";
    //int sOrder = 1;
    //int tOrder = 1;
    int serial_ref_levels = 0;
    int parallel_ref_levels = 0;
    //bool visualization = true;
    double dt = 1.0e-12;
    double dtsf = 0.95;
    //double ti = 0.0;
    //double ts = 1.0;
    //double tf = 40.0;
 
    Array<int> abcs;
    Array<int> dbcs;
 
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                    "Mesh file to use.");
//   args.AddOption(&sOrder, "-so", "--spatial-order",
//                  "Finite element order (polynomial degree).");
//   args.AddOption(&tOrder, "-to", "--temporal-order",
//                  "Time integration order.");
   args.AddOption(&serial_ref_levels, "-rs", "--serial-ref-levels",
                  "Number of serial refinement levels.");
   args.AddOption(&parallel_ref_levels, "-rp", "--parallel-ref-levels",
                  "Number of parallel refinement levels.");
//   args.AddOption(&dtsf, "-sf", "--dt-safety-factor",
//                  "Used to reduce the time step below the upper bound.");
//   args.AddOption(&ti, "-ti", "--initial-time",
//                  "Beginning of time interval to simulate (ns).");
//   args.AddOption(&tf, "-tf", "--final-time",
//                  "End of time interval to simulate (ns).");
//   args.AddOption(&ts, "-ts", "--snapshot-time",
//                  "Time between snapshots (ns).");
//   args.AddOption(&ds_params_, "-ds", "--dielectric-sphere-params",
//                  "Center, Radius, and Permittivity of Dielectric Sphere");
//   args.AddOption(&ms_params_, "-ms", "--magnetic-shell-params",
//                  "Center, Inner Radius, Outer Radius, and Permeability "
//                  "of Magnetic Shell");
//   args.AddOption(&cs_params_, "-cs", "--conductive-sphere-params",
//                  "Center, Radius, and Conductivity of Conductive Sphere");
//   args.AddOption(&dp_params_, "-dp", "--dipole-pulse-params",
//                  "Axis End Points, Radius, Amplitude, "
//                  "Pulse Center (ns), Pulse Width (ns)");
    args.AddOption(&abcs, "-abcs", "--absorbing-bc-surf",
                    "Absorbing Boundary Condition Surfaces");
    args.AddOption(&dbcs, "-dbcs", "--dirichlet-bc-surf",
                    "Dirichlet Boundary Condition Surfaces");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                    "--no-visualization",
                    "Enable or disable GLVis visualization.");
    args.Parse();
    if (!args.Good())
    {
        if (Mpi::Root())
        {
            args.PrintUsage(cout);
        }
        return 1;
    }
    if (Mpi::Root())
    {
        args.PrintOptions(cout);
    }
 
    // Read the (serial) mesh from the given mesh file on all processors.  We can
    // handle triangular, quadrilateral, tetrahedral, hexahedral, surface and
    // volume meshes with the same code.
    Mesh *mesh;
    ifstream imesh(mesh_file);
    if (!imesh)
    {
        if (Mpi::Root())
        {
            cerr << "\nCan not open mesh file: " << mesh_file << '\n' << endl;
        }
        return 2;
    }
    mesh = new Mesh(imesh, 1, 1);
    imesh.close();
    
    // Refine the serial mesh on all processors to increase the resolution. In
    // this example we do 'ref_levels' of uniform refinement.
    for (int l = 0; l < serial_ref_levels; l++)
    {
        mesh->UniformRefinement();
    }
 
    // Define a parallel mesh by a partitioning of the serial mesh. Refine this
    // mesh further in parallel to increase the resolution. Once the parallel
    // mesh is defined, the serial mesh can be deleted.
    ParMesh pmesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
 
    // Refine this mesh in parallel to increase the resolution.
    for (int l = 0; l < parallel_ref_levels; l++)
    {
        pmesh.UniformRefinement();
    }


    // program
    

    return 0;
}
