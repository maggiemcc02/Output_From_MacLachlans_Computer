

# In this file I have my solution discovery code
# This code makes use of all initial guesses in my initial guess pool
# to discover solutions

# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



# imports

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
from coupling import *




def solution_discovery(mesh, uni_grid, guess_list, monitor_func, \
                                                           eps, N, damping, alpha, power, found_sols, found_mesh, call_num, boor_tol, physical_tol, u_0, u_n): 
  
  
  results = []


  global solutions
  global sol_mesh




  for m in range(len(guess_list)):


    solutions = []
    sol_mesh = []


    for i in range(len(found_sols)):

      solutions.append(found_sols[i])
      sol_mesh.append(found_mesh[i])




    # damping = True

    print()
    print('_'*100)
    print('_'*100)
    print('TRYING INITIAL GUESS ' + str(m+1))
    print('_'*100)
    print('_'*100)
    print()

    guess = guess_list[m]


    attempt = ultimate_coupling_unigrid(mesh, uni_grid, guess, monitor_func, eps, N, damping, alpha, power, solutions, sol_mesh, boor_tol, physical_tol, u_0, u_n)

    results.append([len(attempt[0]), attempt[0], attempt[1]]) # number of solutions, solutions, meshes



    # update our found solution list


    
  max_sols = 0
  max_ind = 100

  print()
  print()

  for i in range(len(results)):

    try_max = results[i][0]


    print('For initial guess', i+1, 'we found', try_max, 'solutions')
    print()

    if try_max > max_sols:

      max_sols = try_max

      max_ind = i

  print('We choose the solution set for initial guess', max_ind + 1)
  print()

  chosen_solutions = results[max_ind][1]
  chosen_meshes = results[max_ind][2]

  for i in range(len(chosen_solutions)):

    plt.plot(chosen_meshes[i], chosen_solutions[i], 'blue')
    plt.plot(chosen_meshes[i], [0 for j in range(len(chosen_meshes[i]))], marker = '|', color = 'darkblue')
    plt.xlabel('mesh')
    plt.ylabel('solution approx')
    plt.title('Solution '+ str(i+1) + ' at eps = '+  str(eps))
    plt.savefig('/output_plots/Solution_' + str(i+1) + " at eps = " + str(eps) + ".pdf")


  # return that chosen solution set


  return chosen_solutions, chosen_meshes



