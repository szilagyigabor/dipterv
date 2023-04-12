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
    // width -> y coord, length -> x coord
    const char *new_mesh_file = "discontinuity.mesh";
    double fw = 4;
    double sw = 2;
    double offset = 0.0;
    double fl = 10;
    double sl = 10;
    int nw = 3;
    int nl = 3;
    int order = 1;
    
    // process command line parameters
    OptionsParser args(argc, argv);
    args.AddOption(&new_mesh_file, "-m", "--mesh-out-file",
                   "Output Mesh file to write.");
    args.AddOption(&nw, "-nw", "--num-elements-width",
                   "Number of elements along width. (Must be at least 3)");
    args.AddOption(&nl, "-nl", "--num-elements-length",
                   "Number of elements along length. (Must be at least 2)");
    args.Parse();
    if (!args.Good())
    {
       args.PrintUsage(cout);
       return 1;
    }
    args.PrintOptions(cout);
    if(nw<3 || nl<2)
    {
        printf("\nError: Number of widthwise elements (-nw)must be at least 3 and mnumber of lengthwise elements (-nl) must be at least 2.\n");
        return 1;
    }
    double wider_width = max(fw,sw), narrower_width = min(fw,sw);
    if(offset >= (wider_width-narrower_width)/2.0)
    {
        printf("\nError: offset absolute value is too large\n");
        return 1;
    }

    
    // set finite element space to quads
    Element::Type el_type = Element::QUADRILATERAL;
    // generate initial mesh
    //Mesh mesh = Mesh::MakeCartesian2D(nw, nl, el_type, 1, max(fw,sw), fl+sl);
    //Mesh mesh = Mesh(2,4,1);
    
    Mesh *mesh = new Mesh(2,0,0);
    
    // calculate node coords and number of nodes
    Array<double> xcoord(nw+1);
    Array<double> ycoord(nl+1);
    bool first_is_wider = fw>sw;
    if(fw < sw) offset *= -1.0;
    double bot_width = (wider_width-narrower_width)/2.0+offset;
    double top_width = (wider_width-narrower_width)/2.0-offset;
    // number of widthwise elements that stick out on the bottom of the wider piece (at least one)
    int num_w_el_bot = max(1, (int)(round((double)nw*bot_width/wider_width)));
    // number of widthwise elements for the narrower piece (at least one)
    int num_w_el_mid = max(1, min(nw-2, (int)(round((double)nw*narrower_width/wider_width))));
    // number of widthwise elements that stick out on the top of the wider piece (at least one)
    int num_w_el_top = nw-num_w_el_bot-num_w_el_mid;
    double dy_bot = bot_width/(double)(num_w_el_bot);
    double dy_mid = narrower_width/(double)(num_w_el_mid);
    double dy_top = top_width/(double)(num_w_el_top);
    // add node y coords that define the bottom part
    double y = -dy_bot;
    for(int i=0; i<num_w_el_bot+1; i++)
    {
        y += dy_bot;
        ycoord[i]=y;
    }
    // add node y coords that define the middle part
    for(int i=num_w_el_bot+1; i<num_w_el_bot+1+num_w_el_mid; i++)
    {
        y += dy_mid;
        ycoord[i]=y;
    }
    // add node y coords that define the top part
    for(int i=num_w_el_bot+1+num_w_el_mid; i<num_w_el_bot+1+num_w_el_mid+num_w_el_top; i++)
    {
        y += dy_top;
        ycoord[i]=y;
    }
    // number of lengthwise elements of the first part
    //int num_l_el_first = max(1, nl*(int)(round(fl/(fl+sl))));
    int num_l_el_first = max(1, (int)(round((double)nl*fl/(fl+sl))));
    //printf("num_l_el_first = %d\n",num_l_el_first);
	// number of lengthwise elements of the second part
    int num_l_el_second = nl-num_l_el_first;
    double dx_first = fl/(double)(num_l_el_first);
    double dx_second = sl/(double)(num_l_el_second);
    // add node x coords that define the first part
    double x = -dx_first;
    for(int i=0; i<num_l_el_first+1; i++)
    {
        x += dx_first;
        xcoord[i]=x;
    }
    // add node x coords that define the second part
    for(int i=num_l_el_first+1; i<num_l_el_first+1+num_l_el_second; i++)
    {
        x += dx_second;
        xcoord[i]=x;
    }

	//printf("%f, %f, %f, %f, %f, %f\n", dy_bot, dy_mid, dy_top, dx_first, dx_second, 1.0/0.0);
	// add middle part for both segments
	Array2D<int> vi(nw+1,nl+1); // vertex indices
	// add bottom edge vertices of the middle section
	for(int col=0; col<nl+1; col++)
	{
		vi[num_w_el_bot][col] =
    		mesh->AddVertex(xcoord[col], ycoord[num_w_el_bot]);
	}
	for(int row=num_w_el_bot; row<num_w_el_bot+num_w_el_mid; row++)
	{
		// add top left vertex of the first element of the row
    	vi[row+1][0] =
			mesh->AddVertex(xcoord[0], ycoord[row+1]);
		// add the rest of the vertices and the elements of the row
		for(int col=0; col<nl; col++)
		{
    		vi[row+1][col+1] =
        		mesh->AddVertex(xcoord[col+1], ycoord[row+1]);
			mesh->AddQuad(vi[row][col], vi[row+1][col],
				vi[row+1][col+1], vi[row][col+1]);
		}
		// add boundary edges on the two ends of the middle section (in the x/lentgh direction)
		// boundary has attribute=2 on the x=0 side (port 1), attr.=3 on the other end (port 2)
		// and attr.=1 on the walls
		//mesh->AddBdrSegment(vi[row][0], vi[row+1][0], 2);
		//mesh->AddBdrSegment(vi[row][nl], vi[row+1][nl], 3);
	}
	vi.Print();
	printf("xcoord:\n");
	xcoord.Print(out, 10);
	printf("ycoord:\n");
	ycoord.Print(out, 10);
	

//    // add vertices for "first" transmission line piece
//    for(int i=0; i<NY; i++)
//    {
//        for(int j=0; j<NX; j++)
//        {
//            vert_ind[i][j] = mesh->AddVertex(i*dy, j*dx, attr);
//        }
//    }
//    // add Quad elements
//    for(int i=0; i<NY-1; i++)
//    {
//        for(int j=0; j<NX-1; j++)
//        {
//            mesh->AddQuad(vert_ind[i][j], vert_ind[i+1][j],
//                    vert_ind[i+1][j+1], vert_ind[i][j+1], attr);
//        }
//    }
//    // add boundary edges on top and bottom
//    for(int i=0; i<NX-1; i++)
//    {
//        mesh->AddBdrSegment(vert_ind[0][i+1], vert_ind[0][i], 2);
//        mesh->AddBdrSegment(vert_ind[NY-1][NX-2-i], vert_ind[NY-1][NX-1-i], 3);
//    }
//    // add boundary edges on right and left side
//    for(int i=0; i<NY-1; i++)
//    {
//        mesh->AddBdrSegment(vert_ind[NY-2-i][0], vert_ind[NY-1-i][0], 4);
//        mesh->AddBdrSegment(vert_ind[i+1][NX-1], vert_ind[i][NX-1], 5);
//    }
//    
    mesh->Finalize();

    ofstream ofs("discontinuity.mesh");
    ofs.precision(8);
    mesh->Print(ofs);
    ofs.close();

    delete mesh;

    return 0;
}
