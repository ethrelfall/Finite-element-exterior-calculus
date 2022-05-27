#This is intended to solve L u = f where L is the Hodge Laplacian and u and f are 1-forms, in the square 2D domain (0,2)^2, with u \cdot n and curl(u) equal to zero on the boundary.
#The analytic solution is u=(0.5*((x-1)^2-1), 0)

from firedrake import *

mesh=SquareMesh(25,25,2,2, quadrilateral=False) #uniform mesh so example self-contained but similar results are found on unstructured meshes

#the choice of solution spaces RT1 and CG1 is one of two structure-preserving cases presented in Finite Element Exterior Calculus by D.N.Arnold, p.89
SRT=FunctionSpace(mesh,"RT",1) #the problem is this works but the output looks grim, and it only converges for polynomial order equal to one
#SRT=VectorFunctionSpace(mesh, "CG",1) #solution looks better if this space is used - why?
SP=FunctionSpace(mesh,"CG",1) 
V = SRT*SP

u, sigma = TrialFunctions(V) 
v, tau = TestFunctions(V)

#this mixed formulation is taken from (4.32) in Finite Element Exterior Calculus by D.N.Arnold
a=( inner(sigma, tau)- inner(u, grad(tau)) +inner(grad(sigma),v)+inner(curl(u), curl(v)))*dx

x, y = SpatialCoordinate(mesh)

f=Function(SRT)
f.interpolate(as_vector([-1,0]))
L=inner(f,v) *dx

g=Function(V)

solve(a == L, g, solver_parameters={'ksp_type': 'gmres', 'pc_type': 'none', 'ksp_gmres_restart': 100, 'ksp_monitor_true_residual': None, 'ksp_view': None, 'ksp_converged_reason': None, 'ksp_rtol': 5e-3})

u, sigma = g.split()
File("FEEC_test_1_u.pvd").write(u)
#File("FEEC_test_1_sigma.pvd").write(sigma)
