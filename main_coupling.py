

# In here I will import everything, set parameters, and run codes
# This will be my main code if you will :)




# import tools

import math
import numpy as np
from scipy.linalg import solve_banded
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# import files


from automatic_solver import *
from classic_arclength import *
from coupling import *
from deflated_continuation_code import *
from deflation import *
from equidistribute import *
from solution_discovery_code import *
from initial_guess_pool import *
from L2_optimal import *
from physical_solver import *
from solution_and_jacobian import *



# set parameters


u_0 = 3/2
u_n = 1/2

E = 0.05
N = 100

power = 2.5
alpha = 1

physical_tol = 1e-8
boor_tol = 1e-8

damping = True

eps_0 = 0.05
eps_f = 0.01
delta_eps = 0.01

# set grids


grid = np.linspace(0, 1, N)
mesh = np.array([grid]).T
uni_grid = np.array([grid]).T


# set guesses


guess_list = []

guesses = [initial_guess_1, initial_guess_2, initial_guess_3, initial_guess_4, initial_guess_5]

for j in range(5):

  initial_guess = guesses[j]


  guess = []
  for i in range(len(mesh)):
    if i != 0:
      if i != (len(mesh)-1):
        val = [ initial_guess(mesh[i][0]) ]
        guess.append(val)
  guess = np.array( guess )


  guess_list.append(guess)



# call solution discovery



chosen_sols, chosen_mesh = solution_discovery(mesh, uni_grid, guess_list, M_calc_optimal_L2, \
                                                           E, N, damping, alpha, power, [], [], call_num = 0)










