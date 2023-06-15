

# In this file I have my continued deflation code



# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# import files

from automatic_solver import *
from deflation import *
from physical_solver import *
from equidistribute import *
from system_and_jacobian import *
from classic_arclength import *
from optimal_arc import *
from L2_optimal import *
from numerical_integration import *
from initial_guess_pool import *
from solution_discovery_code import *
from coupling import *






def continued_deflation(eps_0, eps_f, delta_eps,  new_sol_eps, new_mesh_eps, guess, grid, uni_grid, monitor_func, boor_tol, physical_tol, N, alpha, power):




  # list to track number of solutions

  results = []


  eps = eps_0 - delta_eps

  call_num = 1

  while eps >= eps_f:

    call_num += 1


    print()
    print()
    print()
    print()
    print('EPSILON = ', eps)
    print()
    print()
    print()
    print()

    global solutions
    global sol_mesh

    solutions = []
    sol_mesh = []


    # clear out the lists 

    sol_eps = new_sol_eps
    mesh_eps = new_mesh_eps

    # set the initial guesses

    guesses = []

    for j in range(len(sol_eps)):


      guess = []
      for l in range(len(sol_eps[0])):
        if l != 0:
          if l != (len(sol_eps[0])-1):

            val = [ sol_eps[j][l][0] ]

            guess.append(val)

      guess = np.array( guess )

      guesses.append(guess)

    
    new_sol_eps = []
    new_mesh_eps = []

    for i in range(len(sol_eps)):

      print()
      print()
      print('-_'*100)
      print('FINDING SOLUTION NUMBER : ', i, 'at eps = ', eps)
      print('_-'*100)
      print()
      print()


      # find a new solution

      new_sol, new_mesh = automatic_solver_interpU(F, JF, mesh_eps[i], uni_grid,\
                                                  guesses[i], monitor_func, solutions,\
                                                  sol_mesh, boor_tol, physical_tol,\
                                                   eps, N, alpha, power, u_0, u_n)
    
      

      if type(new_sol) == str:

        print('We could not find a solution so we move on')
        print()
        continue
      

      # save to list 

      new_sol_eps.append(new_sol)
      new_mesh_eps.append(new_mesh)
      
      # plot the solution

      # plt.plot(new_mesh, new_sol, 'red')
      # plt.plot(new_mesh, [0 for i in range(len(new_mesh))], marker = "|", color = 'darkred')
      # plt.xlabel('Grid Points')
      # plt.ylabel('Solution Function Approximation')
      # plt.title("Solution Number " + str(i+1) + " at eps = " +str(eps))
      # plt.show()



    # deflate the solutions out and try and find more

    print()
    print()
    print()
    print('NOW GOING TO DEFLATE')
    print()
    print()
    print()



    # solutions = new_sol_eps
    # sol_mesh = new_mesh_eps



    new_sol_deflate, new_mesh_deflate = solution_discovery(mesh, uni_grid, guess_list, monitor_func, \
                                                           eps, N, damping, alpha, power, new_sol_eps, new_mesh_eps, call_num, boor_tol, physical_tol, u_0, u_n)
    


     
    
    # if len(new_sol_deflate) != 0:


      # for def_ind in range(len(new_sol_deflate)):

        # print()
        # print()
    
        # plt.plot(new_mesh_deflate[def_ind], new_sol_deflate[def_ind], 'blue')
        # plt.plot(new_mesh_deflate[def_ind], [0 for i in range(len(new_mesh_deflate[def_ind]))], marker = "|", color = 'darkblue')
        # plt.xlabel('Grid Points')
        # plt.ylabel('Solution Function Approximation')
        # plt.title('Solution '+str(def_ind + 1)+' Found With Deflation')
        # plt.show()


        # print()
        # print()
        # print()


        # # check the solution 


        # mesh1 = []
        # sol1 = []
        # for i in range(len(new_mesh_deflate[def_ind])):
        #   if i != 0: # if not first entry
        #     if i != (len(new_mesh_deflate[def_ind])-1): # if not last entry
        #       mesh1.append([ new_mesh_deflate[def_ind][i][0] ]) # then add it to new vector
        #       sol1.append([new_sol_deflate[def_ind][i][0]])
        # mesh1 = np.array(mesh1)
          

        # # Calculate the mesh width (h) for each grid point
        # # h_i = x_i - x_i-1

        # h_list1 = []
        # for i in range(1, len(new_mesh_deflate[def_ind])):
        #     h = new_mesh_deflate[def_ind][i][0] - new_mesh_deflate[def_ind][i-1][0]
        #     h_list1.append( [h] )
        # mesh_space_1 = np.array(h_list1)


        # check1 = F(sol1, mesh1, mesh_space_1, k, u_0, u_n, eps )


        # if norm(check1) > physical_tol:

        #   print()
        #   print('*'*100)
        #   print('THIS SOLUTION IS NOT GOOD???')
        #   print('*'*100)
        #   print()
        #   print()


        # print('CHECKING THE SOLUTION :', norm(check1))
        # print()
        # print()

    print()
    print()
    print("WE FOUND", len(new_sol_eps), 'SOLUTIONS')
    print()
    print()

    results.append([eps, len(new_sol_eps), new_sol_eps, new_mesh_eps])

    # update epsilon

    eps -= delta_eps


  # print the numbers of solutions

  print()
  print()
  print('_'*100)
  print()

  for i in range(len(results)):

    print('For eps = ', eps, 'we found', results[i][0], 'solutions')
    print()

  print()
  print('_'*100)
  print()
  print()

  

  return results
