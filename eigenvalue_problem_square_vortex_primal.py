from firedrake import *
from firedrake.petsc import PETSc
try:
	from slepc4py import SLEPc
except ImportError:
		import sys
		warning("Unable to import SLEPc, eigenvalue computation not possible (try firedrake-update --slepc)")
		sys.exit(0)

# attempts solve for harmonics of \nabla^2 on an annular domain, primal CG
# fails (either cohomology not reproduced or problems with re-entrant corner divergences)
# finds spurious \lambda \neq 0 lowest eigenvalue - true lowest eigenvalue is zero

mesh=Mesh("eigenvalue_problem_square_vortex.msh")

order=1

V = VectorFunctionSpace(mesh, "CG", order)

eigenmodes_real=Function(V)
eigenmodes_imag=Function(V)

u = TrialFunction(V)
v = TestFunction(V)

x,y = SpatialCoordinate(mesh)

a=(curl(u)*curl(v)+div(u)*div(v))*dx
m = inner(u,v)*dx

# specify essential BCs u \cdot n = 0 (u \times n = 0 is natural)
bc1 = DirichletBC(V.sub(0), 0.0, {21,23,25,27})  # edges with horizontal normals
bc2 = DirichletBC(V.sub(1), 0.0, {20,22,24,26})  # edges with vertical normals

petsc_a = assemble(a, bcs=[bc1,bc2]).M.handle
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

cnt = V.dof_dset.size
print("num dofs = ")
print(cnt)
npa2 = npa1[0:cnt]
npa3 = npa1[0:cnt]

SS = FunctionSpace(mesh, "CG", order)
eigenmodex = Function(SS)
eigenmodey = Function(SS)
npa4 = npa1[::2]
npa5 = npa1[1::2]
eigenmodex.vector()[:] = npa4  # (think this is x-cpt) 
eigenmodey.vector()[:] = npa5  # (""            y""  )

File("eigenvalue_problem_square_vortex_primal_u.pvd").write(eigenmodex, eigenmodey)

print("eigenvalue value:")
print(lam)
