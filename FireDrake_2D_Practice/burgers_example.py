# In this file we will follow the FireDrake tutorial to solve Burgers Equation
# Burgers equation is a nonlinear equation for the advection and diffusion of momentum
# It is given as
# u_t + (u dot grad)u - scalar grad^2u = 0 on the unit square with (n dot grad)u = 0 on the boundary of the unit square



# import firedrake

from firedrake import *

# set the mesh

n = 30
mesh = UnitSquareMesh(n,n)

# Set the function space V
# we choose degree 2 continuous Lagrange polynomials

V = VectorFunctionSpace(mesh, "CG", 2)

# we need a piecewise linear space for output purposes

V_out = VectorFunctionSpace(mesh, "CG", 1)

# since this is a nonlinear problem we don't define trial functions (Newton vibes so we arent filling some set u vector)
# Rather, the trail functions will be created automatically when the residual is differentiated by the nonlinear solver

# we need solution functions for the current and next time step

u_ = Function(V, name = "Velocity")
u = Function(V, name = "VelocityNext")

# define the test function

v = TestFunction(V)


# set the initial condition

x = SpatialCoordinate(mesh)
ic = project(as_vector([sin(pi*x[0]), 0]), V) # project vector of sin values onto function space V

# we start with the current value of u set to the initial condition
# We also use the initial condition as our starting guess for the next value of u

u_.assign(ic)
u.assign(ic)

# set nu (scalar)

nu = 0.0001

# set the timestep

timestep = 1.0/n


# define the residual of the equation

# In the advection term we need to contract the test function v with (u dot grad)u which is the derivative of velocity in the direction u.
# This directional derivative can be written as dot(u, nabla_grad(u)) since nabla_grad(u)[i,j] = dell_i u_j


F = (inner((u - u_)/timestep,v) + inner(dot(u, nabla_grad(u)), v) + nu*inner(grad(u), grad(v))) * dx

# outfile


outfile = File("burgers.pvd")

outfile.write(project(u, V_out, name = 'Velocity')) # needed because output only supports visualization of linear fields so we project onto a linear space


# loop over the time steps
# solve the equation each time and output each result
# firedrakes default used - applying full LU as a preconditioner 

t = 0.0

end = 0.5

while (t <= end):

    solve(F==0, u)

    u_.assign(u)

    t += timestep

    

    try:
        import matplotlib.pyplot as plt
    except:
        warning("Matplotlib not imported")

    try:

        fig, axes = plt.subplots()
        print('First line ok')
        print('plt ok')
        colours = tripcolor(project(u, V_out, name = 'Velocity'), axes = axes)
        print('contours ok')
        print('axes ok')
        print('Second line ok')
        fig.colorbar(colours)
        fig.suptitle("solution Approximation")
        print('Third line ok')
        print()

    except Exception as e:
        warning('Cannot plot figure. Error msg "%s"' %e)

    try:
        plt.show()

    except Exception as e:
        warning("Cannot show figure. Error msg: '%s'" %e)




    outfile.write(project(u, V_out, name = 'Velocity'))

