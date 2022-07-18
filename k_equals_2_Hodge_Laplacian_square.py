from firedrake import *	

# solves k=2 (i.e. 2-form) Hodge Laplace problem on unit square
# this corresponds to \nabla^2 u = f with homogeneous Dirichlet BCs
# note that here sigma is a 1-form and u is a "scalar" (i.e. single-cpt, actually 2-form in 2D),
# i.e. opposite way round to k=1 Hodge Laplacian examples

# this is based on "Finite Element Exterior Calculus: from Hodge Theory to Numerical Stability"
# (arXiv:0906.4325v3) by Arnold, Falk, Winther (AFW); S.2.3.1 p.14 on, outputs Fig.2.2

mesh = UnitSquareMesh(20,20, diagonal='right')

# note re what works (notation means vector func space, scalar func space choices):
# RT1, DG0 is AFW's "working" plot (Fig.2.2 RHS)
# CG1, DG0 is AFW's "fail" plot (Fig.2.2 LHS)
# CGN, CGN gives nonsense
# RTN, CGN ""
# CGN, CG(N-1) seems to work OK
# RTN, CG(N-1) ""

SV = FunctionSpace(mesh, "RT", 1)  # RT not RTE, order should be one more than that of scalar space
#SV=VectorFunctionSpace(mesh, "CG", 1)
SS = FunctionSpace(mesh, "DG", 0)  # zero-order allowed for DG, not CG (it's obvious why as CG0 would imply globally constant)
#SS = FunctionSpace(mesh, "CG", 1)

V = SV * SS

sigma, u = TrialFunctions(V)
tau, v = TestFunctions(V)

a = ( inner(sigma, tau) - u*div(tau) + div(sigma)*v )*dx

f=Function(SS)

x,y = SpatialCoordinate(mesh)
f.interpolate(2*y*(1-y)+2*x*(1-x))  # this was reverse-engineered from AFW's analytic solution u = x(1-x)y(1-y)

L=( f*v )*dx

g=Function(V)
solve(a==L, g)

sigma, u = g.split()  # opposite way round to k=1 Hodge Laplace examples
File("k_equals_2_Hodge_Laplacian_square.pvd").write(u)  # u or sigma
