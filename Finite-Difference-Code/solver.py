'''solver.py
根据所选格式计算数值解，计算l2误差，并绘制动画。计算格式精度时，自动生成多组网格，返回计算的精度p。
'''


import numpy as np
from scheme import FTCS, Lax, Lax_Wendroff, upwind, upwind2, Warming_Beam, leapfrog, Adam_Bashforth
from drawing import animation, plot_solution, plot_l2error
from meshgrid import gridupdate

# 主要计算单元，调用对应格式，返回数值解
def calculate(scheme, u0, x, t, dx, dt, nt, tmax, c, plot=True):
    step = nt - 1
    u_exact = exact(x, dt, step)

    if scheme == 'FTCS':
        u = FTCS(u0, c, step)
    elif scheme == 'Lax':
        u = Lax(u0, c, step)
    elif scheme == 'Lax_Wendroff':
        u = Lax_Wendroff(u0, c, step)
    elif scheme == 'upwind':
        u = upwind(u0, c, step)
    elif scheme == 'upwind2':
        u = upwind2(u0, c, step)
    elif scheme == 'Warming_Beam':
        u = Warming_Beam(u0, c, step)
    elif scheme == 'leapfrog':
        u = leapfrog(u0, c, step)
    elif scheme == 'Adam_Bashforth':
        u = Adam_Bashforth(u0, c, step)
    else:
        raise TypeError("The scheme is not in list")
    if plot:
        animation(x, t, u, dt, c, scheme, u_exact)
        plot_solution(x, t, u, nt, c, scheme, u_exact)

    return np.array(u), np.array(u_exact)

# 理论精度
def theory(scheme):
    if scheme == 'FTCS':
        return None
    elif scheme == 'Lax':
        return 1
    elif scheme == 'Lax_Wendroff':
        return 2
    elif scheme == 'upwind':
        return 1
    elif scheme == 'upwind2':
        return None
    elif scheme == 'Warming_Beam':
        return 2
    elif scheme == 'leapfrog':
        return 2
    elif scheme == 'Adam_Bashforth':
        return 2
    else:
        raise TypeError("The scheme is not in list")


# 生成精确解
def exact(x, dt, step):
    U = []
    t = 1
    while t <= step:
        U.append(np.sin(2 * np.pi * (x - dt * 1 * t)))
        t += 1
    return U

# 计算不同网格尺度的u，比较格式精度
def ordercheck(scheme, nhalf, check, a, xmin, xmax, nx, tmax, nt):
    nx_list, nt_list, grid_list, tmax_list = gridupdate(nhalf, xmin, xmax, nx, tmax, nt)
    l2error_list = []
    dx_list = []
    dt_list = []
    for i in range(0, nhalf):
        nx = nx_list[i]
        nt = nt_list[i]
        tmax = tmax_list[i]
        x, t, dx, dt = grid_list[i]
        u0 = np.sin(2 * np.pi * x)
        c = a * dt / dx
        u, u_exact = calculate(scheme, u0, x, t, dx, dt, nt, tmax, c, plot=False)
        l2err = l2error(u, u_exact, dx)
        l2error_list.append(l2err)
        dx_list.append(dx)
        dt_list.append(dt)
        print(f"c ={c:.3f}, l2err={l2err:.3f}, when nx={nx}, dx={dx}, nt={nt}, dt={dt:.5f}")
    p = p_cal(l2error_list, nhalf)
    plot_l2error(dx_list, l2error_list, scheme, p, check, Theory=theory(scheme))

    return


# 计算l2误差
def l2error(u, ue, dx):
    err = u[0] - ue[0]
    l2err = np.sqrt(np.sum(err ** 2 * dx))
    return l2err

# 计算p,格式精度
def p_cal(l2error_list, nhalf):
    p_list = []
    for i in range(0, nhalf - 1):
        p = np.log(l2error_list[i] / l2error_list[i + 1]) / np.log(2)
        p_list.append(p)
    return np.mean(np.array(p_list))
