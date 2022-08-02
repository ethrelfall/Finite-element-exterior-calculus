from firedrake import *

# solves \nabla^2 u = f, where f is (1,0)
# note this corresponds to the example in Fig.5.1 of the textbook Finite Element Exterior Calculus, D.N. Arnold

mesh=Mesh("re-entrant_corner.msh")

# there are two "good" choices of vector finite element space:

# 1) Raviart-Thomas edge elements - order should be same as that of space SS below
SV=FunctionSpace(mesh,"RTE",3)  # note element type N1curl is equivalent to RTE

# 2) Brezzi-Douglas-Marini edge elements - order should be one less than that of space SS below
#SV=FunctionSpace(mesh, "BDME",2)  # note element type N2curl is equivalent to BDME

# ... there is also the "obvious" choice of CG function space, which does not work well:
#SV=VectorFunctionSpace(mesh, "CG", 3)

SS=FunctionSpace(mesh,"CG",3)  # order should be 1, 2, 3, ...
V = SV*SS

# mixed formulation with sigma = div(u) 
u, sigma = TrialFunctions(V)
v, tau   = TestFunctions(V)

a=( inner(sigma, tau) - inner(u, grad(tau)) +inner(grad(sigma),v)+inner(curl(u),curl(v)) ) *dx

x, y = SpatialCoordinate(mesh)

f=Function(SV)
f.interpolate(as_vector([-1,0]))
L=inner(f,v) *dx

g=Function(V)

solve(a == L, g)

u, sigma = g.split()
File("re-entrant_corner_u.pvd").write(u)
#File("re-entrant_corner_FEEC_sigma.pvd").write(sigma)  # if sigma desired as output

# plot rhs from evaluation of Laplacian of the solution
r = Function(SV)
r.interpolate(grad(div(u))-curl(curl(u)))
File("re-entrant_corner_rhs.pvd").write(r)
