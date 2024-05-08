# Here we do mixed possion tutorial

from firedrake import *
import matplotlib.pyplot as plt

# set mesh

N = 20
mesh = UnitSquareMesh(N,N)


# choose LBB-stable finite element pair for W = V x Sigma

rt = FiniteElement("Raviart-Thomas", triangle, 2, variant = "integral")
Sigma = FunctionSpace(mesh, rt)
V = FunctionSpace(mesh, "DG", 1)


# build the mixed space

W = Sigma * V

sigma, u = TrialFunctions(W)
tau, v = TestFunctions(W)


# set up to solve

x,y = SpatialCoordinate(mesh)

f = 10*exp(-100 * ((x-0.5)**2 + (y-0.5)**2))


a = dot(sigma, tau)*dx + div(tau) * u * dx + div(sigma)*v * dx

L = -f * v * dx


# solve

wh = Function(W) # hold solution

solve( a == L, wh, solver_parameters = {"ksp_type":"minres", "pc_type":"none"})

#plots

from firedrake.pyplot import quiver, tripcolor

sigmah, uh = wh.subfunctions

fig,axes = plt.subplots(ncols = 2, sharex = True, sharey = True)

quiver(sigmah, axes = axes[0])
axes[0].set_aspect('equal')
axes[0].set_title("$\Sigma$")

tripcolor(uh, axes = axes[1])
axes[1].set_aspect("equal")
axes[1].set_title("$u$")

plt.tight_layout()
plt.show()
