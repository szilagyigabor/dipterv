// Mesh for a 2D model of a transmission line discontinuity.
// The transmission line pieces have a rectangular cross-section,
// the first piece has width "fl" and length "fw", ("sl" and "sw" for the secon one).

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    // default values for command line parameters
    const char *new_mesh_file = "discontinuity.mesh";
    double fw = 0.04;
    double sw = 0.02;
    double fl = 0.1;
    double sl = 0.1;
    int nw = 8;
    int nl = 30;
    int order = 1;
    
    // process command line parameters
    OptionsParser args(argc, argv);
    args.AddOption(&new_mesh_file, "-m", "--mesh-out-file",
                   "Output Mesh file to write.");
    args.AddOption(&nw, "-nw", "--num-elements-width",
                   "Number of elements along width.");
    args.AddOption(&nl, "-nl", "--num-elements-length",
                   "Number of elements along length.");
    args.AddOption(&order, "-o", "--mesh-order",
                   "Order (polynomial degree) of the mesh elements.");
    args.Parse();
    if (!args.Good())
    {
       args.PrintUsage(cout);
       return 1;
    }
    args.PrintOptions(cout);
    
    // set finite element space to quads
    Element::Type el_type = Element::QUADRILATERAL;
    // generate initial mesh
    Mesh mesh = Mesh::MakeCartesian2D(nw, nl, el_type, 1, max(fw,sw), fl+sl);
    //Mesh mesh = Mesh(2,4,1);
    mesh.Finalize();
    

    ofstream ofs(new_mesh_file);
    ofs.precision(8);
    mesh.Print(ofs);
    ofs.close();

    return 0;
}
