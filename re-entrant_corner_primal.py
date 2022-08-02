from firedrake import *

# attempt to solve \nabla^2 u = f, where f is (1,0),
# primal method fails (because candidate func space cannot represent singularity)
# and converges to the wrong answer
# see re-entrant_corner.py for a correct implementation
# note this corresponds to the example in Fig.5.1 of the textbook Finite Element Exterior Calculus, D.N. Arnold

mesh=Mesh("re-entrant_corner_primal.msh")

V=VectorFunctionSpace(mesh,"CG",1)

# primal formulation
u=TrialFunction(V)
v=TestFunction(V)

a=(curl(u)*curl(v)+div(u)*div(v))*dx

f=Function(V)
x, y = SpatialCoordinate(mesh)
f.interpolate(as_vector([1,0]))

L=inner(f,v)*dx

g=Function(V)

# specify essential BCs u \cdot n = 0 (u \times n = 0 is natural)
bc1 = DirichletBC(V.sub(0), 0.0, {15,17,19})  # edges with horizontal normals
bc2 = DirichletBC(V.sub(1), 0.0, {16,18,20})  # edges with vertical normals

solve(a==L,g, bcs=[bc1, bc2])

File("re-entrant_corner_primal.pvd").write(g)

# plot the RHS from evaluating Laplacian of the solution - surprisingly, the equation seems well-satisfied away from the re-entrant corner (NOTE use at least order-2 for this)
r = Function(V)
r.interpolate(grad(div(g))-curl(curl(g)))
File("re-entrant_corner_primal_rhs.pvd").write(r)
