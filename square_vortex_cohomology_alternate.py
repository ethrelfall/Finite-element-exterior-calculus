from firedrake import *

# an alternate method for finding the cohomology of the square domain with non-centred square hole
# same solution as zero-eigenvalue mode in eigenvalue_problem_square_vortex.py
# problem is solved here without using FEEC techniques
# method: solve for the difference between the usual (-y/r^2, x/r^2) vortex and the square case,
# the difference is the gradient of a scalar which satisfies a Laplacian Neumann problem
# cf. deforming a contour by incorporating a topologically-trivial loop
# note: if comparing output with that of eigenvalue_problem_square_vortex.py, be aware normalization likely to differ

mesh=Mesh("square_vortex_cohomology_alternate.msh")

SS=FunctionSpace(mesh,"CG",3)

u = TrialFunction(SS)
v = TestFunction(SS)

x,y = SpatialCoordinate(mesh)

a = inner(grad(u), grad(v)) *dx  # scalar Laplacian

# Neumann BC
Fxp=Function(SS) # xp = largest x-coord domain edge
Fxp.interpolate(-y/(x**2+y**2))
Fxm=Function(SS) # xm = smallest x-coord domain edge
Fxm.interpolate( y/(x**2+y**2))
Fyp=Function(SS)
Fyp.interpolate( x/(x**2+y**2))
Fym=Function(SS)
Fym.interpolate(-x/(x**2+y**2))

L = inner(Fym,v)*ds(20)+inner(Fxp,v)*ds(21)+inner(Fyp,v)*ds(22)+inner(Fxm,v)*ds(23)+inner(Fyp,v)*ds(24)+inner(Fxm,v)*ds(25)+inner(Fym,v)*ds(26)+inner(Fxp,v)*ds(27)

g=Function(SS)

nullspace = VectorSpaceBasis(constant=True)

solve(a == L, g, nullspace=nullspace)

# now calculate the "flow field" by subtracting the circular vortex from the gradient of the scalar
# the subtraction makes the boundary conditions right for the square domain
SV=VectorFunctionSpace(mesh, "CG", 2)
flow = Function(SV)
flow.interpolate(as_vector([y/(x**2+y**2),-x/(x**2+y**2)])+grad(g))

File("square_vortex_cohomology_alternate_flow.pvd").write(flow)
