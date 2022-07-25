# vector Helmholtz equation on square with boundary conditions

from firedrake import *

mesh=SquareMesh(25,25,2,2, quadrilateral=False)

SV = FunctionSpace(mesh,"RTE",1)
SS = FunctionSpace(mesh,"CG",1) 
V = SV*SS

u, sigma = TrialFunctions(V) 
v, tau = TestFunctions(V)

a=( inner(sigma, tau)- inner(u, grad(tau)) +inner(grad(sigma),v)+inner(curl(u), curl(v))+inner(u,v))*dx  # is Helmholtz


f=Function(SV)
f.interpolate(as_vector([-1,1]))
L=inner(f,v) *dx

g=Function(V)

bc = DirichletBC(V.sub(0), as_vector([0.0,0.0]), "on_boundary")  # this constrains both cpts not just tangential
solve(a==L,g,bcs=bc)

u, sigma = g.split()
File("square_bc_test_submit_u.pvd").write(u)

