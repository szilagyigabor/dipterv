#include "pod_solver.hpp"
#include <stdio.h>

using namespace std;
using namespace mfem;
//using namespace mfem::common;

int main(int argc, char *argv[])
{
    
    Mpi::Init();
    Hypre::Init();
 
    const char *mesh_file = "discontinuity.mesh";
    int sOrder = 1;
    //int tOrder = 1;
    int serial_ref_levels = 0;
    int parallel_ref_levels = 0;
    bool visualization = true;
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
    args.AddOption(&sOrder, "-so", "--spatial-order",
                   "Finite element order (polynomial degree).");
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
    Mesh *serial_mesh;
    ifstream imesh(mesh_file);
    if (!imesh)
    {
        if (Mpi::Root())
        {
            cerr << "\nCan not open mesh file: " << mesh_file << '\n' << endl;
        }
        return 2;
    }
    serial_mesh = new Mesh(imesh, 1, 1);
    imesh.close();
    
    // Refine the serial mesh on all processors to increase the resolution. In
    // this example we do 'ref_levels' of uniform refinement.
    for (int l = 0; l < serial_ref_levels; l++)
    {
        serial_mesh->UniformRefinement();
    }
 
    // Define a parallel mesh by a partitioning of the serial mesh. Refine this
    // mesh further in parallel to increase the resolution. Once the parallel
    // mesh is defined, the serial mesh can be deleted.
    ParMesh mesh(MPI_COMM_WORLD, *serial_mesh);
    serial_mesh->Clear();
 
    // Refine this mesh in parallel to increase the resolution.
    for (int l = 0; l < parallel_ref_levels; l++)
    {
        mesh.UniformRefinement();
    }

    /********************** program *********************/
    
    // FE space for magnetix field
    ND_FECollection B_fec(sOrder, mesh.Dimension());
    ParFiniteElementSpace B_fespace(&mesh, &B_fec);
    HYPRE_BigInt B_num_dofs = B_fespace.GlobalTrueVSize();

    // FE space for magnetic field
    L2_FECollection E_fec(sOrder, mesh.Dimension());
    ParFiniteElementSpace E_fespace(&mesh, &E_fec);
    HYPRE_BigInt E_num_dofs = E_fespace.GlobalTrueVSize();

    if (Mpi::Root())
    {
        cout << "Number of E field unknowns: " << E_num_dofs << endl
            << "Number of B field unknowns: " << B_num_dofs << endl;
    }
    
    // get E field boundaries
    Array<int> E_dirichlet_boundary_dofs, E_neumann_1_boundary_dofs, E_neumann_2_boundary_dofs;
    E_fespace.GetBoundaryTrueDofs(E_dirichlet_boundary_dofs, 1); // 1: dirichlet boundaries
    E_fespace.GetBoundaryTrueDofs(E_neumann_1_boundary_dofs, 2); // 2: first port
    E_fespace.GetBoundaryTrueDofs(E_neumann_2_boundary_dofs, 3); // 3: second port
    
    // get B field boundaries
    // maybe not needed, for now, no absrobing boundaries
    Array<int> B_dirichlet_boundary_dofs, B_neumann_1_boundary_dofs, B_neumann_2_boundary_dofs;
    B_fespace.GetBoundaryTrueDofs(B_dirichlet_boundary_dofs, 1); // 1: dirichlet boundaries
    B_fespace.GetBoundaryTrueDofs(B_neumann_1_boundary_dofs, 2); // 2: first port
    B_fespace.GetBoundaryTrueDofs(B_neumann_2_boundary_dofs, 3); // 3: second port

    // 

    /******************* program vége *******************/

    return 0;
}
