



from firedrake import *
import numpy as np
import matplotlib.pyplot as plt


# import the mesh

mesh0 = Mesh('my_polygon.msh')



def Mesh_Iter(u_sol, initial_mesh, V, iter):



  # set the initial_mesh

  mesh = initial_mesh 





  for i in range(iter):




    # Compute the mesh density function

    grad_u = grad(u_sol)
    u_x = grad_u[0]
    u_y = grad_u[1]
    M = sqrt(1 + u_x**2 + u_y**2)


    # indicate a new mesh solve


    print()
    print()
    print('~' * 35)
    print('MESH SOLVE ', i+1)
    print('~'*35)
    print()
    print()



    # # solve the mesh problem  - x




    # set the test and and trial functions

    x  = TrialFunction(V)
    v = TestFunction(V)

    # define the  problem

    a =  ( M *inner( grad(x), grad(v) ) )* dx

    L = ( Constant(0) * v  ) * dx

    # define the BC's

    bc1 =  DirichletBC(V, Constant(0), 1)
    bc2 = DirichletBC(V, Constant(1), 2)
    bc3 = DirichletBC(V, xi, 3)
    bc4 = DirichletBC(V, xi, 4)

    bcs = [ bc1, bc2, bc3 , bc4]

    # solve

    x_sol = Function(V)

    solve(a == L, x_sol, bcs = bcs)




    # solve the mesh problem - y




    # set the test and and trial functions

    y  = TrialFunction(V)
    v = TestFunction(V)

    # define the nonlinear problem

    a = (M *inner( grad(y), grad(v) ) )* dx

    L = ( Constant(0) * v  ) * dx

    # define the BC's

    bc1 =  DirichletBC(V, eta, 1)
    bc2 = DirichletBC(V, eta, 2)
    bc3 = DirichletBC(V, Constant(0), 3)
    bc4 = DirichletBC(V, Constant(1), 4)

    bcs = [ bc1, bc2, bc3, bc4]

    # solve

    y_sol = Function(V)

    solve(a == L, y_sol, bcs = bcs)



    # Create and plot the new mesh

    N = np.shape(x_sol.dat.data)[0]
    mesh_values = np.zeros([N,2])
    mesh_values[:,0] = x_sol.dat.data
    mesh_values[:,1] = y_sol.dat.data
    mesh.coordinates.dat.data[:] = mesh_values
    fig, axes = plt.subplots()
    triplot(mesh, axes=axes)
    plt.title("Adapted Mesh for the Current Solution")
    plt.show()

    print()
    print()


    # interpolate the solution onto a function space defined by the new mesh


    V2 = FunctionSpace(mesh, "CG", 1)

    u_interp = interpolate(u_sol, V2)

    fig, axes = plt.subplots()
    levels = np.linspace(0, 1, 20)
    contours_y = tricontour(u_sol,  axes = axes, levels = levels, colors = 'red', label = 'last solution')
    #contours_x = tricontour(u_interp,  axes = axes, levels = levels, colors = 'lightpink', label = "interpolated solution")
    fig.colorbar(contours_y)
    #fig.colorbar(contours_x)
    #plt.legend()
    fig.suptitle("NonInterpolated Solution")

    fig, axes = plt.subplots()
    levels = np.linspace(0, 1, 20)
    #contours_y = tricontour(u_sol,  axes = axes, levels = levels, colors = 'red', label = 'last solution')
    contours_x = tricontour(u_interp,  axes = axes, levels = levels, colors = 'green', label = "interpolated solution")
    #fig.colorbar(contours_y)
    fig.colorbar(contours_x)
    #plt.legend()
    fig.suptitle("Interpolated Solution")

    fig, axes = plt.subplots()
    levels = np.linspace(0, 1, 20)
    contours_y = tricontour(u_sol,  axes = axes, levels = levels, colors = 'red', label = 'last solution')
    contours_x = tricontour(u_interp,  axes = axes, levels = levels, colors = 'green', label = "interpolated solution")
    fig.colorbar(contours_y)
    fig.colorbar(contours_x)
    #plt.legend()
    fig.suptitle("Both For Comparison")


    u_sol = u_interp

    V = V2



  

  return mesh











# SOLVE THE PHYSICAL PROBLEM


#mesh0 = UnitSquareMesh(20, 20)
sigma = Constant(0.25)


# set the function space - a vector function space

V = FunctionSpace(mesh0, "CG", 1)

# Access the Mesh Coordinates

xi, eta = SpatialCoordinate(mesh0) # mesh is \xi and \eta

# set the test and and trial functions

u  = TrialFunction(V)
v = TestFunction(V)

# set f

f = Function(V)

gauss = exp((-((xi - 1/2)**2 + (eta - 1/2)**2)) / (2 * (sigma)**2))

f.interpolate( (-1/(sigma**2)) * gauss * ( 2 - (1/(sigma**2))*(xi - 1/2)**2  - (1/(sigma**2))*(eta - 1/2)**2 ))

# define the problem

a = ( -1 * inner( grad(u), grad(v) ) ) * dx # NEGATIVE

L = ( inner( f, v ) ) * dx

# define the BC's

bcy1 = exp(-((0 - 1/2)**2 + (eta - 1/2)**2) / (2*sigma**2)) # x = 0
bcy2 = exp(-((1 - 1/2)**2 + (eta - 1/2)**2) / (2*sigma**2)) # x = 1
bcx1 = exp(-((xi - 1/2)**2 + (0 - 1/2)**2) / (2*sigma**2)) # y = 0
bcx2 = exp(-((xi - 1/2)**2 + (1 - 1/2)**2) / (2*sigma**2)) # y = 1

bc1 =  DirichletBC(V, bcy1 , 1)
bc2 = DirichletBC(V, bcy2 , 2)
bc3 = DirichletBC(V, bcx1 , 3)
bc4 = DirichletBC(V, bcx2 , 4)

bcs = [ bc1, bc2, bc3, bc4]

# solve the problem

u_sol = Function(V)

solve(a == L, u_sol, bcs = bcs)

# Colour Plot

fig, axes = plt.subplots()
colors = tripcolor(u_sol, axes = axes)
plt.title('Solution Approximation')
fig.colorbar(colors)
plt.show()
print()
print()




# compute the exact solution

u_exact = Function(V).interpolate( exp((-((xi - 1/2)**2 + (eta - 1/2)**2)) / (2 * (sigma)**2)) )

# compute the error

error = errornorm(u_exact, u_sol)

print("Error is: ",  error)
print()





# now mesh solve

new_mesh = Mesh_Iter(u_sol, mesh0, V, 10)






