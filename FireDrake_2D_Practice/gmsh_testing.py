# In this file we do the tutorial for generating meshes
# In particular, we do a problem with used linear Lagrangain finite elements


# import firedrake

from firedrake import *

#  load the mesh

mesh = Mesh('immersed_domain_msh')

# define the function space

V = FunctionSpace(mesh, "CG", 1)

# define the trial and test functions

u = TrialFunction(V)
v = TestFunction(V)


# define the bilinear form of the problem

a = 2*dot(grad(v), grad(u))*dx(4) + dot(grad(v), grad(u))*dx(3) + v*u*dx

# the tag 4 is the disc tag and 3 is the rectangle without the disc

# define the linear form

L = Constant(5.) * v * dx + Constant(3.) * v('+') * dS(13)

# set the homogeneous Dirichlet bc's on the rectangle boundaries
# tag 11 - horizontal edges, tag 12 - vertical edges

DirBC = DirichletBC(V, 0, [11, 12])

# define u to hold the solution

u = Function(V)


# solve

solve(a==L, u, bcs = DirBC, solver_parameters = {'ksp_type':'cg'})
