'''meshgrid.py
划分网络，计算精度时返回多套网格。
'''

import numpy as np


# 生成网格
def grid(xmin, xmax, nx, tmax, nt):
    dx = (xmax - xmin) / (nx - 1)
    dt = tmax / (nt - 1)
    x = np.linspace(xmin, xmax, nx)
    t = np.linspace(0, tmax, nt)
    return x, t, dx, dt


def gridupdate(nhalf, xmin, xmax, onx, tmax, ont):
    nx_list = []
    nt_list = []
    grid_list = []
    tmax_list = []

    nx = onx
    for i in range(0, nhalf):
        nx_list.append(nx)
        nx = int(nx * 2 - 1)

    nt = ont
    for i in range(0, nhalf):
        nt_list.append(nt)
        nt = int(nt * 2 - 1)

    for i in range(0, nhalf):
        x, t, dx, dt = grid(xmin, xmax, nx_list[i], tmax, nt_list[i])
        grid_list.append([x, t, dx, dt])
        tmax_list.append(tmax)
        # tmax = tmax * 2

    return nx_list, nt_list, grid_list, tmax_list
