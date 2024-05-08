# in here I solve poissons equation with dirichlet conditions


# import firedrake, create mesh, set function space

from firedrake import *
import matplotlib.pyplot as plt
from firedrake.pyplot import tripcolor


mesh = UnitSquareMesh(10,10)
V = FunctionSpace(mesh, "Lagrange", 1)


# set exact solution and rhs of pde

x,y = SpatialCoordinate(mesh)
u_exact = sin(pi*x)*sin(pi*y)
f = 2 * pi ** 2 * u_exact


# test and trial functions

u = TrialFunction(V)
v = TestFunction(V)

# bilinear and linea

a = dot(grad(u), grad(v)) * dx
L = f * v * dx

# set the zero Dirichlet conditions on the edge of the domain

boundary_ids = (1,2,3,4)
bcs = DirichletBC(V, 0, boundary_ids)


# solve

uh = Function(V)
solve(a==L, uh, bcs = bcs, solver_parameters = {'ksp_type':'cg', 'pc_type':'none'})


# plot 

fig,axes = plt.subplots()
collection = tripcolor(uh, axes = axes)
fig.colorbar(collection)
plt.show()




