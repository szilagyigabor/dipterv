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
    // width -> y coord, length -> x coord, height -> z coord
    const char *new_mesh_file = "discontinuity.mesh";
    double fw = 4.0;
    double sw = 2.0;
    double os = 0.0;
    double fl = 7.0;
    double sl = 7.0;
    double h = 3.0;
    int nw = 3;
    int nl = 2;
    int nh = 1;
    
    // process command line parameters
    OptionsParser args(argc, argv);
    args.AddOption(&new_mesh_file, "-m", "--mesh-out-file",
                   "Output Mesh file to write.");
    args.AddOption(&nw, "-nw", "--num-elements-width",
                   "Number of elements along width. (Must be at least 3)");
    args.AddOption(&nl, "-nl", "--num-elements-length",
                   "Number of elements along length. (Must be at least 2)");
    args.AddOption(&nh, "-nh", "--num-elements-height",
                   "Number of elements along height.");
    args.AddOption(&h, "-h", "--height",
                   "Height of model.");
    args.AddOption(&fw, "-fw", "--first-width",
                   "Width of first section.");
    args.AddOption(&sw, "-sw", "--second-width",
                   "Width of second section.");
    args.AddOption(&os, "-os", "--offset",
                   "Offset of the second element from the center.");
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
    if(os >= (max(fw,sw)-min(fw,sw))/2.0)
    {
        printf("\nError: os absolute value is too large\n");
        return 1;
    }
    GeomParams gp(fw, sw, os, fl, sl, h, nw, nl, nh);
    
    // set finite element space to quads
    Element::Type el_type = Element::HEXAHEDRON;
    // generate initial mesh
    Mesh *mesh = new Mesh(3,0,0);
    
	// add middle part for both segments
	Array2D<int> vi_a(nw+1,nl+1); // vertex indices for layer "a"
	Array2D<int> vi_b(nw+1,nl+1); // vertex indices for layer "b"
	
    // add bottom vertices of the middle section
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
		mesh->AddBdrSegment(vi[row][0], vi[row+1][0], 2);
		mesh->AddBdrSegment(vi[row+1][nl], vi[row][nl], 3);
	}

	// add the top and bottom sections to the wider part
	// determine where is the wider part
	int wide_start, wide_end;
	int narrow_start, narrow_end;
	int attr_start, attr_end;
	if(first_is_wider)
	{
		wide_start = 0;
		wide_end = num_l_el_first;
		narrow_start = num_l_el_first;
		narrow_end = nl;
		attr_start = 2;
		attr_end = 1;
	}
	else
	{
		wide_start = num_l_el_first;
		wide_end = nl;
		narrow_start = 0;
		narrow_end = num_l_el_first;
		attr_start = 1;
		attr_end = 3;
	}
    // add elements on top of the wider part
	for(int row=num_w_el_bot+num_w_el_mid; row<nw; row++)
	{
		// add top left vertex of the first element of the row
    	vi[row+1][0] =
			mesh->AddVertex(xcoord[0], ycoord[row+1]);
		// add the rest of the vertices and the elements of the row
		for(int col=wide_start; col<wide_end; col++)
		{
    		vi[row+1][col+1] =
        		mesh->AddVertex(xcoord[col+1], ycoord[row+1]);
			mesh->AddQuad(vi[row][col], vi[row+1][col],
				vi[row+1][col+1], vi[row][col+1]);
		}
		// add boundary edges on the two ends of the top section (in the x/lentgh direction)
		// boundary has attribute=2 on the x=0 side (port 1), attr.=3 on the other end (port 2)
		// and attr.=1 on the walls
		mesh->AddBdrSegment(vi[row][wide_start], vi[row+1][wide_start], attr_start);
		mesh->AddBdrSegment(vi[row][wide_end], vi[row+1][wide_end], attr_end);
	}
	// add elements on bottom of the wider part
	for(int row=num_w_el_bot; row>0; row--)
	{
		// add bottom left vertex of the first element of the row
    	vi[row-1][0] =
			mesh->AddVertex(xcoord[0], ycoord[row-1]);
		// add the rest of the vertices and the elements of the row
		for(int col=wide_start; col<wide_end; col++)
		{
    		vi[row-1][col+1] =
        		mesh->AddVertex(xcoord[col+1], ycoord[row-1]);
			mesh->AddQuad(vi[row][col], vi[row-1][col],
				vi[row-1][col+1], vi[row][col+1]);
		}
		// add boundary edges on the two ends of the top section (in the x/lentgh direction)
		// boundary has attribute=2 on the x=0 side (port 1), attr.=3 on the other end (port 2)
		// and attr.=1 on the walls
		mesh->AddBdrSegment(vi[row][wide_start], vi[row-1][wide_start], attr_start);
		mesh->AddBdrSegment(vi[row][wide_end], vi[row-1][wide_end], attr_end);
	}

	// add boundary edges of the top and bottom segments
	for(int col=wide_start; col<wide_end; col++)
	{
		mesh->AddBdrSegment(vi[nw][col], vi[nw][col+1], 1);
		mesh->AddBdrSegment(vi[0][col], vi[0][col+1], 1);
	}

	// add boundaries on the top and bottom edge of the narrower part
	for(int col=narrow_start; col<narrow_end; col++)
	{
		mesh->AddBdrSegment(vi[num_w_el_bot][col],
			vi[num_w_el_bot][col+1], 1);
		mesh->AddBdrSegment(vi[num_w_el_bot+num_w_el_mid][col],
			vi[num_w_el_bot+num_w_el_mid][col+1], 1);
	}

    mesh->FinalizeMesh();

    ofstream ofs("discontinuity.mesh");
    ofs.precision(8);
    mesh->Print(ofs);
    ofs.close();

    delete mesh;

    return 0;
}

class GeomParams
{
public:
    double wider_width;//ok
    double narrower_width;//ok
    int num_w_el_bot;//ok
    int num_w_el_mid;//ok
    int num_w_el_top;//ok
    int num_l_el_first;//ok
    int num_l_el_second;//ok
    int num_h_el;
    int wide_start;//ok
    int wide_end;//ok
    int narrow_start;//ok
    int narrow_end;//ok
    int attr_start;//ok
    int attr_end;//ok
    double bot_width; //ok
    double top_width; //ok
    bool first_is_wider; //ok

    GeomParams(double, double, double, double, double, double, int, int, int);
    double GetX(int idx){ return xcoords[idx]; }
    double GetY(int idx){ return ycoords[idx]; }
    double GetZ(int idx){ return zcoords[idx]; }

private:
    Array<int> xcoords;
    Array<int> ycoords;
    Array<int> zcoords;
}

void AddElementLayer(Mesh *mesh, Vector2D &vi_lower, Vector2D &vi_upper, GeomParams *gp, int z_idx)
{
    
}

void AddVertexLayer(Mesh *mesh, Vector2D &vi, GeomParams *gp, int z_idx)
{
    int nl = gp->num_l_el_first + gp->num_l_el_second;
    // add bottom vertices of the middle section
	for(int col=0; col<nl+1; col++)
	{
		vi[gp->num_w_el_bot][col] =
    		mesh->AddVertex(gp->GetX(col), gp->GetY(gp->num_w_el_bot);
	}
	for(int row=gp->num_w_el_bot; row<gp->num_w_el_bot+gp->num_w_el_mid; row++)
	{
		// add top left vertex of the first element of the row
    	vi[row+1][0] =
			mesh->AddVertex(gp->GetX(0), gp->GetY(row+1);
		// add the rest of the vertices and the elements of the row
		for(int col=0; col<nl; col++)
		{
    		vi[row+1][col+1] =
        		mesh->AddVertex(gp->GetX(col+1), gp->GetY(row+1));
			mesh->AddQuad(vi[row][col], vi[row+1][col],
				vi[row+1][col+1], vi[row][col+1]);
		}
		// add boundary edges on the two ends of the middle section (in the x/lentgh direction)
		// boundary has attribute=2 on the x=0 side (port 1), attr.=3 on the other end (port 2)
		// and attr.=1 on the walls
		//mesh->AddBdrSegment(vi[row][0], vi[row+1][0], 2);
		//mesh->AddBdrSegment(vi[row+1][nl], vi[row][nl], 3);
        //TODO
	}

	// add the top and bottom sections to the wider part
	// determine where is the wider part
	int wide_start, wide_end;
	int narrow_start, narrow_end;
	int attr_start, attr_end;
	if(first_is_wider)
	{
		wide_start = 0;
		wide_end = num_l_el_first;
		narrow_start = num_l_el_first;
		narrow_end = nl;
		attr_start = 2;
		attr_end = 1;
	}
	else
	{
		wide_start = num_l_el_first;
		wide_end = nl;
		narrow_start = 0;
		narrow_end = num_l_el_first;
		attr_start = 1;
		attr_end = 3;
	}
    // add elements on top of the wider part
	for(int row=num_w_el_bot+num_w_el_mid; row<nw; row++)
	{
		// add top left vertex of the first element of the row
    	vi[row+1][0] =
			mesh->AddVertex(xcoord[0], ycoord[row+1]);
		// add the rest of the vertices and the elements of the row
		for(int col=wide_start; col<wide_end; col++)
		{
    		vi[row+1][col+1] =
        		mesh->AddVertex(xcoord[col+1], ycoord[row+1]);
			mesh->AddQuad(vi[row][col], vi[row+1][col],
				vi[row+1][col+1], vi[row][col+1]);
		}
		// add boundary edges on the two ends of the top section (in the x/lentgh direction)
		// boundary has attribute=2 on the x=0 side (port 1), attr.=3 on the other end (port 2)
		// and attr.=1 on the walls
		mesh->AddBdrSegment(vi[row][wide_start], vi[row+1][wide_start], attr_start);
		mesh->AddBdrSegment(vi[row][wide_end], vi[row+1][wide_end], attr_end);
	}
	// add elements on bottom of the wider part
	for(int row=num_w_el_bot; row>0; row--)
	{
		// add bottom left vertex of the first element of the row
    	vi[row-1][0] =
			mesh->AddVertex(xcoord[0], ycoord[row-1]);
		// add the rest of the vertices and the elements of the row
		for(int col=wide_start; col<wide_end; col++)
		{
    		vi[row-1][col+1] =
        		mesh->AddVertex(xcoord[col+1], ycoord[row-1]);
			mesh->AddQuad(vi[row][col], vi[row-1][col],
				vi[row-1][col+1], vi[row][col+1]);
		}
		// add boundary edges on the two ends of the top section (in the x/lentgh direction)
		// boundary has attribute=2 on the x=0 side (port 1), attr.=3 on the other end (port 2)
		// and attr.=1 on the walls
		mesh->AddBdrSegment(vi[row][wide_start], vi[row-1][wide_start], attr_start);
		mesh->AddBdrSegment(vi[row][wide_end], vi[row-1][wide_end], attr_end);
	}

	// add boundary edges of the top and bottom segments
	for(int col=wide_start; col<wide_end; col++)
	{
		mesh->AddBdrSegment(vi[nw][col], vi[nw][col+1], 1);
		mesh->AddBdrSegment(vi[0][col], vi[0][col+1], 1);
	}

	// add boundaries on the top and bottom edge of the narrower part
	for(int col=narrow_start; col<narrow_end; col++)
	{
		mesh->AddBdrSegment(vi[num_w_el_bot][col],
			vi[num_w_el_bot][col+1], 1);
		mesh->AddBdrSegment(vi[num_w_el_bot+num_w_el_mid][col],
			vi[num_w_el_bot+num_w_el_mid][col+1], 1);
	}
}

GeomParams::GeomParams(double fw, double sw, double os, double fl, double sl, double h, int nw, int nl, int nh)
    : num_h_el(nh)
{
    wider_width = max(fw,sw);
    narrower_width = min(fw,sw);
    first_is_wider = fw>sw;
    if(!first_is_wider) os *= -1.0;
    bot_width = (wider_width-narrower_width)/2.0+os;
    top_width = (wider_width-narrower_width)/2.0-os;
    // number of widthwise elements that stick out on the bottom of the wider piece (at least one)
    num_w_el_bot = max(1, min(nw-2, (int)(round((double)nw*bot_width/wider_width))));
    // number of widthwise elements for the narrower piece (at least one)
    num_w_el_mid = max(1, min(nw-2, (int)(round((double)nw*narrower_width/wider_width))));
    // number of widthwise elements that stick out on the top of the wider piece (at least one)
    num_w_el_top = nw-num_w_el_bot-num_w_el_mid;
    // number of lengthwise elements of the first part
    num_l_el_first = max(1, min(nl-1, (int)(round((double)nl*fl/(fl+sl)))));
	// number of lengthwise elements of the second part
    num_l_el_second = nl-num_l_el_first;
	if(first_is_wider)
	{
		wide_start = 0;
		wide_end = num_l_el_first;
		narrow_start = num_l_el_first;
		narrow_end = nl;
		attr_start = 2;
		attr_end = 1;
	}
	else
	{
		wide_start = num_l_el_first;
		wide_end = nl;
		narrow_start = 0;
		narrow_end = num_l_el_first;
		attr_start = 1;
		attr_end = 3;
	}

    // generate Y coordinates
    double dy_bot = bot_width/(double)(num_w_el_bot);
    double dy_mid = narrower_width/(double)(num_w_el_mid);
    double dy_top = top_width/(double)(num_w_el_top);
    // add node y coords that define the bottom part
    double y = 0.0;
    for(int i=0; i<num_w_el_bot+1; i++)
    {
        ycoord[i]=y;
        y += dy_bot;
    }
    // add node y coords that define the middle part
    y = ycoord[num_w_el_bot];
    for(int i=num_w_el_bot+1; i<num_w_el_bot+1+num_w_el_mid; i++)
    {
        y += dy_mid;
        ycoord[i]=y;
    }
    // add node y coords that define the top part
    y = ycoord[num_w_el_bot+num_w_el_mid];
    for(int i=num_w_el_bot+1+num_w_el_mid; i<num_w_el_bot+1+num_w_el_mid+num_w_el_top; i++)
    {
        y += dy_top;
        ycoord[i]=y;
    }

    // generate X coordinates
    double dx_first = fl/(double)(num_l_el_first);
    double dx_second = sl/(double)(num_l_el_second);
    // add node x coords that define the first part
    double x = 0.0;
    for(int i=0; i<num_l_el_first+1; i++)
    {
        xcoord[i]=x;
        x += dx_first;
    }
    // add node x coords that define the second part
    x = xcoord[num_l_el_first];
    for(int i=num_l_el_first+1; i<num_l_el_first+1+num_l_el_second; i++)
    {
        x += dx_second;
        xcoord[i]=x;
    }

    // generate Z coordinates
    double dz = h/(double)(nh);
    double z = 0.0;
    for(int i=0; i<num_h_el+1; i++)
    {
        z += dz;
        zcoord[i]=z;
    }
}
