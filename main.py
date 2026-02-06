import numpy as np
import matplotlib.pyplot as plt

j_grid = 1000       # number of grid points
n_grid = 3001       # number of time steps n

C = 3/2            # CFL condition = c*Δt/Δy

source_1 = 400      # source location of EM wave Case 1
source_2 = 10       # source location of EM wave Case 2
source = source_1

Ez = np.zeros(j_grid)       # Electric field (E_z)
Hx = np.zeros(j_grid - 1)   # Magnetic field (H_x)

Ez_prev_left = 0.0      # store previous value of Ez on the left
Ez_prev_right = 0.0     # store previous value of Ez on the right

for n in range(n_grid):                     # loop over time
    for j in range(j_grid -1):              # loop over space in defined grid
        Hx[j] -= (Ez[j+1] - Ez[j]) * C      # update Hx based on previous Ez
    
    for j in range(1, j_grid - 1):
        Ez[j] -= (Hx[j] - Hx[j-1]) * C      # update Ez based on previous Hx

    # Boundary condition for absorbing (Mur boundary condition)
    Ez[0] = Ez_prev_left + (C - 1) / (C + 1) * (Ez[1] - Ez[0])      # define for left boundary
    Ez_prev_left = Ez[1]                                            # store for next iteration
    Ez[-1] = Ez_prev_right + (C - 1) / (C + 1) * (Ez[-2] - Ez[-1])  # define for right boudary
    Ez_prev_right = Ez[-2]

    source_function_1 = np.sin(2 * np.pi * n * C / 20)        # source funtion Case 1 
    source_function_2 = np.exp(- ((n / 30) - 5) ** 2)         # source funtion Case 2

    Ez[source] = source_function_1

    if n == 500 or n == 1000 or n== 3000:       # PLOT at timesteps 500, 1000 and 3000
        plt.figure()
        plt.plot(Ez)
        plt.title(f'Timestep {n} CFL = {C}')
        plt.xlabel('Grid point')
        plt.ylabel('Ez')
        if Ez[source] == source_function_1:
            plt.suptitle('Case 1')
            plt.savefig(f'case1_CFL{C}_time{n}.png') 
        elif Ez[source] == source_function_2:
            plt.suptitle('Case 2')
            plt.savefig(f'case2_CFL{C}_time{n}.png') 
plt.show()