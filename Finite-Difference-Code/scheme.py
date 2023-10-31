''' scheme.py
保存了8个有限差分格式和周期边界条件设置
'FTCS','Lax','Lax_Wendroff','upwind','upwind2','Warming_Beam','leapfrog','Adam_Bashforth'
'''
import numpy as np
import copy


# FTCS格式
def FTCS(u0, r, step):
    u = InitPeriodicBoundary(u0, 'FTCS')
    t = 0
    U = []
    while t < step:
        u = PeriodicBoundary(u, 'FTCS')
        u[1:-1] = u[1:-1] - 0.5 * r * (u[2:] - u[:-2])
        U.append(u.copy()[1:-1])
        t += 1
    return U


# Lax格式
def Lax(u0, r, step):
    u = InitPeriodicBoundary(u0, 'Lax')
    t = 0
    U = []
    while t < step:
        u = PeriodicBoundary(u, 'Lax')
        u[1:-1] = 0.5 * (1 - r) * u[2:] + 0.5 * (1 + r) * u[:-2]
        U.append(u.copy()[1:-1])
        t += 1
    return U


# Lax-Wendroff格式
def Lax_Wendroff(u0, r, step):
    u = InitPeriodicBoundary(u0, 'Lax_Wendroff')
    t = 0
    U = []
    while t < step:
        u = PeriodicBoundary(u, 'Lax_Wendroff')
        u[1:-1] = u[1:-1] - 0.5 * r * (u[2:] - u[:-2]) + 0.5 * r ** 2 * (u[2:] - 2 * u[1:-1] + u[:-2])
        U.append(u.copy()[1:-1])
        t += 1
    return U


# 迎风格式
def upwind(u0, r, step):
    u = InitPeriodicBoundary(u0, 'upwind')
    t = 0
    U = []
    while t < step:
        u = PeriodicBoundary(u, 'upwind')
        u[1:] = u[1:] - r * (u[1:] - u[:-1])
        U.append(u.copy()[1:])
        t += 1
    return U


# 空间二阶迎风格式
def upwind2(u0, r, step):
    u = InitPeriodicBoundary(u0, 'upwind2')
    t = 0
    U = []
    while t < step:
        u = PeriodicBoundary(u, 'upwind2')
        u[2:] = u[2:] - 0.5 * r * (3 * u[2:] - 4 * u[1:-1] + u[:-2])
        U.append(u.copy()[2:])
        t += 1
    return U


# 迎风Warming-Beam格式
def Warming_Beam(u0, r, step):
    u = InitPeriodicBoundary(u0, 'Warming_Beam')
    t = 0
    U = []
    while t < step:
        u = PeriodicBoundary(u, 'Warming_Beam')
        u[2:] = u[2:] - r * (u[2:] - u[1:-1]) - 0.5 * r * (1 - r) * (u[2:] - 2 * u[1:-1] + u[:-2])
        U.append(u.copy()[2:])
        t += 1
    return U


# 蛙跳格式
def leapfrog(u0, r, step):
    u = InitPeriodicBoundary(u0, 'leapfrog')
    u_0 = copy.deepcopy(u)
    u_1 = copy.deepcopy(u_0)
    u_1 = PeriodicBoundary(u_1, 'leapfrog')
    u_1[1:-1] = u_1[1:-1] - 0.5 * r * (u_1[2:] - u_1[:-2]) + 0.5 * r ** 2 * (u_1[2:] - 2 * u_1[1:-1] + u_1[:-2])

    t = 0
    Ujump = [copy.deepcopy(u_0), copy.deepcopy(u_1)]
    U = []
    while t < step:
        u_0 = PeriodicBoundary(Ujump[t], 'leapfrog')
        u_1 = PeriodicBoundary(Ujump[t + 1], 'leapfrog')
        u[1:-1] = u_0[1:-1] - r * (u_1[2:] - u_1[:-2])
        Ujump.append(copy.deepcopy(u))
        t += 1
        U.append(u.copy()[1:-1])
    return U


# Adam-Bashforth格式
def Adam_Bashforth(u0, r, step):
    u = InitPeriodicBoundary(u0, 'leapfrog')
    u_0 = copy.deepcopy(u)
    u_1 = copy.deepcopy(u_0)
    u_1 = PeriodicBoundary(u_1, 'leapfrog')
    u_1[1:-1] = u_1[1:-1] - 0.5 * r * (u_1[2:] - u_1[:-2]) + 0.5 * r ** 2 * (u_1[2:] - 2 * u_1[1:-1] + u_1[:-2])

    t = 0
    Ujump = [copy.deepcopy(u_0), copy.deepcopy(u_1)]
    U = []
    while t < step:
        u_0 = PeriodicBoundary(Ujump[t], 'Adam_Bashforth')
        u_1 = PeriodicBoundary(Ujump[t + 1], 'Adam_Bashforth')
        u[1:-1] = u[1:-1] - 0.25 * r * (3 * (u_1[2:] - u_1[:-2]) - (u_0[2:] - u_0[:-2]))
        Ujump.append(copy.deepcopy(u))
        t += 1
        U.append(u.copy()[1:-1])
    return U


def PeriodicBoundary(u, scheme):
    if scheme in ['FTCS', 'Lax', 'Lax_Wendroff', 'leapfrog', 'Adam_Bashforth']:
        ufull = u
        ufull[0] = u[-2]
        ufull[-1] = u[1]
    elif scheme in ['upwind']:
        ufull = u
        ufull[0] = u[-1]
    elif scheme in ['upwind2', 'Warming_Beam']:
        ufull = u
        ufull[0] = u[-2]
        ufull[1] = u[-1]
    else:
        raise TypeError('格式错误')

    return ufull


def InitPeriodicBoundary(u0, scheme):
    if scheme in ['FTCS', 'Lax', 'Lax_Wendroff', 'leapfrog', 'Adam_Bashforth']:
        u = np.zeros(u0.shape[0] + 2)
        u[1:-1] = u0
    elif scheme in ['upwind']:
        u = np.zeros(u0.shape[0] + 1)
        u[1:] = u0
    elif scheme in ['upwind2', 'Warming_Beam']:
        u = np.zeros(u0.shape[0] + 2)
        u[2:] = u0

    else:
        raise TypeError('格式错误')
    return u
