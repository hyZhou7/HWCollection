'''
2023.10.31
采用不同格式求解一节波动方程
u_t+u_x=0
'''

import numpy as np
from meshgrid import grid
from solver import calculate, ordercheck

#  step0: initial setting
a = 1
xmin = 0
xmax = 3
# ---------------------------------setting now----------------------------------
tmax = 4
schemelib = ['FTCS', 'Lax', 'Lax_Wendroff', 'upwind', 'upwind2', 'Warming_Beam', 'leapfrog', 'Adam_Bashforth']

# step1: grid setting
nx = 101
nt = 101

# step2: scheme setting
scheme = schemelib[0]

# ---------------------------------runing now----------------------------------
x, t, dx, dt = grid(xmin, xmax, nx, tmax, nt)
u0 = np.sin(2 * np.pi * x)
c = a * dt / dx
print("c =", c)

# calculate error
calculate(scheme, u0, x, t, dx, dt, nt, tmax, c, plot=True)

# check accuracy
check = 'x'
nhalf = 4
if check is not None:
    print(f"Start check the {scheme} accuracy.")
    ordercheck(scheme, nhalf, check, a, xmin, xmax, nx, tmax, nt)

print("End of script: Solve a section of wave equations using different schemes")

