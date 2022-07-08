from firedrake import *

# an alternate method for finding the cohomology of the square domain with two circular holes
# same solutions as the two zero-eigenvalue modes in eigenvalue_problem_double_vortex.py
# problem is solved here without using FEEC techniques
# method: solve for the difference between counter-rotating double vortex / co-rotating double vortex and this case,
# the difference is the gradient of a scalar which satisfies a Laplacian Neumann problem
# cf. deforming a contour by incorporating a topologically-trivial loop
# note: if comparing output with that of eigenvalue_problem_double_vortex.py, be aware normalization likely to differ,
# and also the eigensolver does not necessarily find the pure counter / co-rotating modes

# note the counter-rotating solution assumes an infinite plane and can be found as grad(\sigma) in bipolar coords,
# note nice Apollonian circles property of streamlines / equipotentials;
# the co-rotating solution used here is a vortex (\Gamma ln(r-r_0)) pair, streamlines non-circular hence
# needs additional Neumann condition on the two inner circles

mesh = Mesh("eigenvalue_problem_double_vortex.msh")

SS = FunctionSpace( mesh, "CG", 3)

u = TrialFunction(SS)
v = TestFunction(SS)

x, y = SpatialCoordinate(mesh)

a = inner(grad(u), grad(v)) *dx  # scalar Laplacian

F31=Function(SS)
F32=Function(SS)
F33=Function(SS)
F34=Function(SS)
F35=Function(SS)
F36=Function(SS)
g=Function(SS)

# 1. counter-rotating case, based on deformation by gradient

# Neumann BC
F31.interpolate( (x**2-y**2-1)/((x**2+y**2)**2+1+2*(y**2-x**2)))  # bot
F32.interpolate( 2*x*y/((x**2+y**2)**2+1+2*(y**2-x**2)))          # RHS
F33.interpolate(-(x**2-y**2-1)/((x**2+y**2)**2+1+2*(y**2-x**2)))  # top
F34.interpolate(-2*x*y/((x**2+y**2)**2+1+2*(y**2-x**2)))          # LHS 

#F35.interpolate(0.0)  # LHS inner circle
#F36.interpolate(0.0)  # RHS inner circle

L = inner(F31,v)*ds(31)+inner(F32,v)*ds(32)+inner(F33,v)*ds(33)+inner(F34,v)*ds(34)+inner(F35,v)*ds(35)+inner(F36,v)*ds(36)

nullspace = VectorSpaceBasis(constant=True)  # Neumann problem

solve(a == L, g, nullspace=nullspace)

# now calculate the "flow field" by subtracting the circular vortex from the gradient of the scalar
# the subtraction makes the boundary conditions right for the domain
SV=VectorFunctionSpace(mesh, "CG", 2)
counterrot_flow = Function(SV)
counterrot_flow.interpolate(as_vector([-2*y*x/((x*x+y*y)*(x*x+y*y)+1+2*(y*y-x*x)),(x*x-y*y-1)/((x*x+y*y)*(x*x+y*y)+1+2*(y*y-x*x))])+grad(g))

File("double_vortex_cohomology_alternate_counterrot.pvd").write(counterrot_flow)

# 2. co-rotating case, based on deformation by gradient

# Neumann BC
F31.interpolate( x*(x**2+y**2-1.0)/((x**2+y**2)**2+1.0+2.0*(y**2-x**2)))
F32.interpolate( y*(x**2+y**2+1.0)/((x**2+y**2)**2+1.0+2.0*(y**2-x**2)))
F33.interpolate(-x*(x**2+y**2-1.0)/((x**2+y**2)**2+1.0+2.0*(y**2-x**2)))
F34.interpolate(-y*(x**2+y**2+1.0)/((x**2+y**2)**2+1.0+2.0*(y**2-x**2)))

F35.interpolate(-y/sqrt((x**2+y**2)**2+1.0+2.0*(y**2-x**2)))
F36.interpolate( y/sqrt((x**2+y**2)**2+1.0+2.0*(y**2-x**2)))

L = inner(F31,v)*ds(31)+inner(F32,v)*ds(32)+inner(F33,v)*ds(33)+inner(F34,v)*ds(34)+inner(F35,v)*ds(35)+inner(F36,v)*ds(36)

nullspace = VectorSpaceBasis(constant=True)  # Neumann problem

solve(a == L, g, nullspace=nullspace)

# now calculate the "flow field" by subtracting the circular vortex from the gradient of the scalar
# the subtraction makes the boundary conditions right for the domain
SV=VectorFunctionSpace(mesh, "CG", 2)
corot_flow = Function(SV)
corot_flow.interpolate(as_vector([-y*(x*x+y*y+1)/((x*x+y*y)*(x*x+y*y)+1+2*(y*y-x*x)),x*(x*x+y*y-1)/((x*x+y*y)*(x*x+y*y)+1+2*(y*y-x*x))])+grad(g))

File("double_vortex_cohomology_alternate_corot.pvd").write(corot_flow)


