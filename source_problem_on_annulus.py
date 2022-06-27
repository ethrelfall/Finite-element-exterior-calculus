from firedrake import *

# solves \nabla^2 u = f, where f is (0,x)
# note this corresponds to the example in Fig.5.2 of the textbook Finite Element Exterior Calculus, D.N. Arnold

mesh=Mesh("source_problem_on_annulus.msh")

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
v, tau = TestFunctions(V)

x, y = SpatialCoordinate(mesh)

a = (inner(sigma,tau)-inner(u, grad(tau))+inner(grad(sigma),v)+inner(curl(u),curl(v))) *dx

f=Function(SV)
f.interpolate(as_vector([0,-x]))

g=Function(V)

# null space: specified to avoid pollution of the output by the cohomology (represented by the usual vortex solution)
# see https://www.firedrakeproject.org/demos/ma-demo.py by C. Cotter for implementation on mixed function space
q=Function(SV)

q.interpolate(as_vector([-y/(x**2+y**2),x/(x**2+y**2)]))  # this could be replaced by a numerical determination of the 1-cohomology
normq = 1/sqrt(assemble(inner(q,q)*dx)) # normalization factor
q*=normq

V_basis = VectorSpaceBasis({q})
nullspace=MixedVectorSpaceBasis(V,[V_basis,V.sub(1)])
transposenullspace = VectorSpaceBasis({})
scalarprod = assemble(inner(f,q)*dx)
f-=q*(scalarprod)  # q normed above
print("normalization check:")
print(assemble(inner(f,q)*dx))  # check orthogonalization - output should be near-zero
L=inner(f,v)*dx

# solver setup: this is copied from C. Cotter's example referenced above (as are the comments)

# We need to set quite a few solver options, so we'll put them into a
# dictionary. ::

sp_it = {

# We'll only use stationary preconditioners in the Schur complement, so
# we can get away with GMRES applied to the whole mixed system ::

#
   "ksp_type": "gmres",

# We set up a Schur preconditioner, which is of type "fieldsplit". We also
# need to tell the preconditioner that we want to eliminate :math:`\sigma`,
# which is field "1", to get an equation for :math:`u`, which is field "0". ::

#
   "pc_type": "fieldsplit",
   "pc_fieldsplit_type": "schur",
   "pc_fieldsplit_0_fields": "1",
   "pc_fieldsplit_1_fields": "0",

# The "selfp" option selects a diagonal approximation of the A00 block. ::

#
   "pc_fieldsplit_schur_precondition": "selfp",

# We just use ILU to approximate the inverse of A00, without a KSP solver, ::

#
   "fieldsplit_0_pc_type": "ilu",
   "fieldsplit_0_ksp_type": "preonly",

# and use GAMG to approximate the inverse of the Schur complement matrix. ::

#
   "fieldsplit_1_ksp_type": "preonly",
   "fieldsplit_1_pc_type": "gamg",
   "fieldsplit_1_mg_levels_pc_type": "sor",

# Finally, we'd like to see some output to check things are working, and
# to limit the KSP solver to 20 iterations. ::

#
   "ksp_monitor": None,
   "ksp_max_it": 200,
   "snes_monitor": None,
   "ksp_rtol": 1.0e-3  # added by ET, tolerance to which solution is solved
   }

solve(a == L, g, nullspace=nullspace, solver_parameters=sp_it)

u, sigma = g.split()
File("source_problem_on_annulus_u.pvd").write(u)
#File("source_problem_on_annulus_sigma.pvd").write(sigma)  # if sigma desired as output
