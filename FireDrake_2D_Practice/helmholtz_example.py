# In this file we will solve the 'Simple Helmholtz Equation' from the FireDrake Tutorials


# problem : -grad^2u + u = f with gradu*n = 0 on boundsry of unit square. We choose f = (1 + 8pi^2)cos(2pix)cos(2piy)

# solution : u = cos(2pix)cos(2piy)




# the 10x10 unit square mesh

from firedrake import *

mesh = UnitSquareMesh(10, 10)

# define the function space

V = FunctionSpace(mesh, "CG", 1) # linear lagrange

# trial functions

u = TrialFunction(V)
v = TestFunction(V)


# declare a function over the function space

f = Function(V)

# Declare f as the RHS 

x, y = SpatialCoordinate(mesh)

f.interpolate((1+8*pi*pi)*cos(x*pi*2)*cos(y*pi*2))


# define the bilinear form of the weak problem

# bilinear : INT(grad(u) inner product grad(v) + u inner product v) # lhs of weak form

a = (inner(grad(u), grad(v)) + inner(u,v)) * dx

# linear : INT( vf )

l = inner(f, v)*dx

# Define u to hold the solution


u = Function(V)


# solve the problem

solve( a == l, u, solver_parameters = {'ksp_type':'cg', 'pc_type': 'none'})

# we instruct PETSc to employ the conjugate gradient method


# save the ouput file

File('helmholtz.pvd').write(u)


# check the L2 norm


#f.interpolate(cos(x*pi*2)*cos(y*pi*2))
#print(sqrt(assemble(dot(u-f, u-f))*dx))



# plot

try:
    import matplotlib.pyplot as plt
except:
    warning("Matplotlib not imported")

try:

    fig, axes = plt.subplots()
    print('First line ok')
    print('plt ok')
    colours = tripcolor(u, axes = axes)
    print('contours ok')
    print('axes ok')
    print('Second line ok')
    fig.colorbar(colours)
    fig.suptitle("solution Approximation of Helmholtz")
    print('Third line ok')
    print()

except Exception as e:
    warning('Cannot plot figure. Error msg "%s"' %e)

try:
    plt.show()

except Exception as e:
    warning("Cannot show figure. Error msg: '%s'" %e)
