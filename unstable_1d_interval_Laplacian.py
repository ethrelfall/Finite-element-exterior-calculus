from firedrake import *

# attempt to solve -u'' = f homog. Dirichlet problem in 1D
# see Fig.1.2 of "Finite element exterior calculus: from Hodge theory to numerical stability"
# Arnold, Falk, Winther, and also their ref.25.

mesh = IntervalMesh(14, -1,1)

SV = VectorFunctionSpace(mesh, "CG", 2) # first-order works, second-order is (expectedly) problematic
SS = FunctionSpace(mesh, "DG", 0)

V = SV*SS

sigma, u = TrialFunctions(V)
tau, v = TestFunctions(V)

a = ( inner(sigma, tau) - u*div(tau) + div(sigma)*v )*dx

f = Function(SS)
x = SpatialCoordinate(mesh)
f.interpolate(0.25*pi*pi*cos(0.5*pi*x[0]))
#f.interpolate(1.0+0.0*x[0]) # a simpler example, see the ref.25 mentioned above

L = (f*v)*dx

g = Function(V)
solve(a == L, g)

sigma, u = g.split()
# note that second-order produces Paraview file that
# crashes Paraview if do plot over line!
# solution is "save data" as csv
# then reconstruct basis functions
File("unstable_1d_interval_Laplacian.pvd").write(sigma) 

