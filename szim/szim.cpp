#include "mfem.hpp"
#include <stdio.h>

#define MU0 1.25663706212e-6
#define EPS0 8.8541878128e-12

using namespace std;
using namespace mfem;
//using namespace mfem::common;

// initial E field based on position "x"
double EFieldFunc(const Vector &);
double JFieldFunc(const Vector &, double);
void BFieldFunc(const Vector &, Vector &);
double muInvFunc(const Vector &);

int main(int argc, char *argv[])
{

	const char *mesh_file = "discontinuity.mesh";
	int order = 1;
	//int tOrder = 1;
	int serial_ref_levels = 0;
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
	args.AddOption(&order, "-o", "--order",
				  "Finite element order (polynomial degree).");
//   args.AddOption(&tOrder, "-to", "--temporal-order",
//				  "Time integration order.");
   args.AddOption(&serial_ref_levels, "-rs", "--serial-ref-levels",
				  "Number of serial refinement levels.");
//   args.AddOption(&dtsf, "-sf", "--dt-safety-factor",
//				  "Used to reduce the time step below the upper bound.");
//   args.AddOption(&ti, "-ti", "--initial-time",
//				  "Beginning of time interval to simulate (ns).");
//   args.AddOption(&tf, "-tf", "--final-time",
//				  "End of time interval to simulate (ns).");
//   args.AddOption(&ts, "-ts", "--snapshot-time",
//				  "Time between snapshots (ns).");
//   args.AddOption(&ds_params_, "-ds", "--dielectric-sphere-params",
//				  "Center, Radius, and Permittivity of Dielectric Sphere");
//   args.AddOption(&ms_params_, "-ms", "--magnetic-shell-params",
//				  "Center, Inner Radius, Outer Radius, and Permeability "
//				  "of Magnetic Shell");
//   args.AddOption(&cs_params_, "-cs", "--conductive-sphere-params",
//				  "Center, Radius, and Conductivity of Conductive Sphere");
//   args.AddOption(&dp_params_, "-dp", "--dipole-pulse-params",
//				  "Axis End Points, Radius, Amplitude, "
//				  "Pulse Center (ns), Pulse Width (ns)");
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
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);
	// Read the (serial) mesh from the given mesh file on all processors.  We can
	// handle triangular, quadrilateral, tetrahedral, hexahedral, surface and
	// volume meshes with the same code.
	ifstream imesh(mesh_file);
	if (!imesh)
	{
		cerr << "\nCan not open mesh file: " << mesh_file << '\n' << endl;
		return 2;
	}
	Mesh mesh(imesh, 1, 1);
	imesh.close();
	mesh.SetCurvature(1); // needed for interpolation at physical points

	// Refine the serial mesh on all processors to increase the resolution. In
	// this example we do 'ref_levels' of uniform refinement.
	for (int l = 0; l < serial_ref_levels; l++)
	{
		mesh.UniformRefinement();
	}


	// program, for now only 2-dimensional meshes
	int dim = 2;

	FunctionCoefficient muInvCoef(muInvFunc); // time-DEPENDENT


	// define E FESpace: H1 continous scalar valued elements
	// also used for J and sigmaE
	H1_FECollection E_fec(order, dim);
	FiniteElementSpace E_fespace(&mesh, &E_fec);

	// define B FESpace: edge elements
	ND_FECollection B_fec(order, dim);
	FiniteElementSpace B_fespace(&mesh, &B_fec);

	// define rotE = (-dB/dt) operator
	MixedBilinearForm rotE(&E_fespace, &B_fespace);
	rotE.AddDomainIntegrator(new MixedScalarWeakCurlIntegrator); // itt lehet egy coeffitient
	rotE.SetAssemblyLevel(AssemblyLevel::FULL);
	rotE.Assemble();
	rotE.Finalize();

	// define rotB = (mu*(J+sigmaE+dD/dt))
	MixedBilinearForm rotB(&B_fespace, &E_fespace);
	rotB.AddDomainIntegrator(new MixedScalarCurlIntegrator(muInvCoef)); // itt is lehet egy coef
	rotB.SetAssemblyLevel(AssemblyLevel::FULL);
	rotB.Assemble();
	rotB.Finalize();

	// grid functions
	GridFunction E(&E_fespace);
	//GridFunction J(&E_fespace);
	GridFunction B(&B_fespace);
	E = 0.0;
	B = 0.0;
//	J = 0.0;

	FunctionCoefficient InitialEFieldCoef(EFieldFunc); // time-independent
//	FunctionCoefficient JFieldCoef(JFieldFunc); // time-DEPENDENT
	VectorFunctionCoefficient InitialBFieldCoef(dim, BFieldFunc); // time-independent

	// set initial fields
	E.ProjectCoefficient(InitialEFieldCoef);
//	J.ProjectCoefficient(JFieldCoef);
	B.ProjectCoefficient(InitialBFieldCoef);

	Vector pos(2);
	pos[0] = 0.1; pos[1] = 0.1;
	cout << "E0(x=0.1,y=0.1)" << EFieldFunc(pos) << endl;
	cout << "E(x=0.1,y=0.1)" << EFieldFunc(pos) << endl;

	FindPointsGSLIB interp;
	interp.Setup(mesh);
	interp.FreeData();
//	printf("%d",E.GetNodes());

	// step time
	Vector e, b;
	E.GetTrueDofs(e);
	B.GetTrueDofs(b);
	Vector dedt(e.Size()), dbdt(b.Size());
	dedt=0.0;
	dbdt=0.0;

	cout << "e.Size(): " << e.Size() << endl;
	cout << "b.Size(): " << b.Size() << endl;
	cout << "rotE.Height(): " << rotE.Height() << endl;
	cout << "rotB.Height(): " << rotB.Height() << endl;
	cout << "rotE.Width(): " << rotE.Width() << endl;
	cout << "rotB.Width(): " << rotB.Width() << endl;

	rotE.Mult(e, dbdt);
	E.SetFromTrueDofs(dedt);
	rotB.Mult(dbdt, dedt);
	B.SetFromTrueDofs(dbdt);

	// save fields
	E.Save("E.gf");
	B.Save("B.gf");

	// save mesh
	mesh.Save("mesh.mesh");

//	delete mesh;
	return 0;
}

// initial E field based on position "x"
double
EFieldFunc(const Vector &x)
{
   return sin(x[0])*sin(x[1]);
}

// J field based on position "x" and time "t"
double
JFieldFunc(const Vector &x, double t)
{
   return 0.0;
}

// initial B field based on position "x"
void
BFieldFunc(const Vector &x, Vector &B)
{
	B.SetSize(2);
	B = 0.0;
//	B[0] = (x[1]-1.0)/10.0;
//	B[1] = (x[0]-1.0)/10.0;
}

// permeability based on position "x"
double
muInvFunc(const Vector &x)
{
	return 1/MU0;
}

