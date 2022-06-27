from firedrake import *
from firedrake.petsc import PETSc
try:
	from slepc4py import SLEPc
except ImportError:
		import sys
		warning("Unable to import SLEPc, eigenvalue computation not possible (try firedrake-update --slepc)")
		sys.exit(0)

# solves for harmonics of \nabla^2 on an annular domain

mesh=Mesh("eigenvalue_problem_single_vortex.msh")

# there are two "good" choices of vector finite element space:

# 1) Raviart-Thomas edge elements - order should be same as that of space SS below
SV=FunctionSpace(mesh,"RTE",3)  # note element type N1curl is equivalent to RTE

# 2) Brezzi-Douglas-Marini edge elements - order should be one less than that of space SS below
#SV=FunctionSpace(mesh, "BDME",2)  # note element type N2curl is equivalent to BDME

SS=FunctionSpace(mesh,"CG",3)  # order should be 1, 2, 3, ...
V = SV*SS

eigenmodes_real=Function(V)
eigenmodes_imag=Function(V)

u, sigma = TrialFunctions(V)
v, tau   = TestFunctions(V)

x,y = SpatialCoordinate(mesh)

a = (inner(sigma, tau)-inner(grad(tau),u)+inner(grad(sigma),v)+inner(curl(u),curl(v)))*dx

m = inner(u,v)*dx
petsc_a = assemble(a).M.handle
petsc_m = assemble(m).M.handle

num_eigenvalues=1

opts = PETSc.Options()
opts.setValue("eps_gen_hermitian", None)
opts.setValue("st_pc_factor_shift_type", "NONZERO")
opts.setValue("eps_type", "krylovschur")
opts.setValue("eps_smallest_real", None)
opts.setValue("eps_tol", 1e-10)

es = SLEPc.EPS().create(comm=COMM_WORLD)
es.setDimensions(num_eigenvalues)
es.setOperators(petsc_a, petsc_m)
es.setFromOptions()

st=es.getST()
st.setType(SLEPc.ST.Type.SINVERT)
es.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)

es.solve()

nconv = es.getConverged()
print("number of converged eigenvalues:")
print(nconv)
vr, vi = petsc_a.getVecs()
lam = es.getEigenpair(0, vr, vi)  #mode 0 spans cohomology space - choose mode here

import numpy as np
npa1 = np.array(vr.getSize())
npa1 = vr

cnt = SV.dof_dset.size
npa2 = npa1[0:cnt]

eigenmode=Function(SV)
eigenmode.vector()[:] = npa2

File("eigenvalue_problem_single_vortex_u.pvd").write(eigenmode)

print("eigenvalue value:")
print(lam)

