// Mesh for a 2D model of a transmission line discontinuity.
// The transmission line pieces have a rectangular cross-section,
// the first piece has width "fl" and length "fw", ("sl" and "sw" for the secon one).

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// number of nodes in x or y directions
#define NX 10
#define NY 12

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
//    args.AddOption(&nw, "-nw", "--num-elements-width",
//                   "Number of elements along width.");
//    args.AddOption(&nl, "-nl", "--num-elements-length",
//                   "Number of elements along length.");
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
    //Mesh mesh = Mesh::MakeCartesian2D(nw, nl, el_type, 1, max(fw,sw), fl+sl);
    //Mesh mesh = Mesh(2,4,1);
    
    Mesh *mesh = new Mesh(2,0,0);
    // rectangle of 2 by 2 elements which have size dx and dy
    // vertices indexed like this:
    // 0,0 0,1 0,2
    // 1,0 1,1 1,2
    // 2,0 2,1 2,2
    double dx = 1.0, dy = 1.0;
    int attr = 1;
    int vert_ind[NY][NX];
    // add vertices
    for(int i=0; i<NY; i++)
    {
        for(int j=0; j<NX; j++)
        {
            vert_ind[i][j] = mesh->AddVertex(i*dy, j*dx, attr);
        }
    }
    // add Quad elements
    for(int i=0; i<NY-1; i++)
    {
        for(int j=0; j<NX-1; j++)
        {
            mesh->AddQuad(
                    vert_ind[i][j], vert_ind[i+1][j], vert_ind[i+1][j+1], vert_ind[i][j+1], attr);
        }
    }
    // add boundary edges on top and bottom
    for(int i=0; i<NX-1; i++)
    {
        mesh->AddBdrSegment(vert_ind[0][i+1], vert_ind[0][i], 2);
        mesh->AddBdrSegment(vert_ind[NY-1][NX-2-i], vert_ind[NY-1][NX-1-i], 3);
    }
    // add boundary edges on right and left side
    for(int i=0; i<NY-1; i++)
    {
        mesh->AddBdrSegment(vert_ind[NY-2-i][0], vert_ind[NY-1-i][0], 4);
        mesh->AddBdrSegment(vert_ind[i+1][NX-1], vert_ind[i][NX-1], 5);
    }
    
    mesh->FinalizeQuadMesh();

    ofstream ofs("discontinuity.mesh");
    ofs.precision(8);
    mesh->Print(ofs);
    ofs.close();

    delete mesh;

    return 0;
}
