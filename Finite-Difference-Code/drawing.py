'''drawing.py
绘制数值解的动画，不同时刻数值解图，格式精度比较图。
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# 生成数值解图
def plot_solution(x, t, u, nt, c, scheme, u_exact):
    fig, ax = plt.subplots(2, 2)

    ax[0][0].plot(x, u[int(nt * 0.05) - 2], color='cornflowerblue', linestyle='-', label="Numer.")
    ax[0][0].plot(x, u_exact[int(nt * 0.05) - 2], color='lightcoral', linestyle='-.', label="Exact")
    ts = t[int(nt * 0.05) - 2]
    ax[0][0].set_title(f"t={ts:.3f}")

    ax[0][1].plot(x, u[int(nt * 0.25) - 2], color='cornflowerblue', linestyle='-', label="Numer.")
    ax[0][1].plot(x, u_exact[int(nt * 0.25) - 2], color='lightcoral', linestyle='-.', label="Exact")
    ts = t[int(nt * 0.25) - 2]
    ax[0][1].set_title(f"t={ts:.3f}")

    ax[1][0].plot(x, u[int(nt * 0.50) - 2], color='cornflowerblue', linestyle='-',label="Numer.")
    ax[1][0].plot(x, u_exact[int(nt * 0.50) - 2], color='lightcoral', linestyle='-.', label="Exact")
    ts = t[int(nt * 0.50) - 2]
    ax[1][0].set_title(f"t={ts:.3f}")

    ax[1][1].plot(x, u[int(nt * 0.75) - 2], color='cornflowerblue', linestyle='-', label="Numer.")
    ax[1][1].plot(x, u_exact[int(nt * 0.75) - 2], color='lightcoral', linestyle='-.', label="Exact")
    ts = t[int(nt * 0.75) - 2]
    ax[1][1].set_title(f"t={ts:.3f}")

    for i, row in enumerate(ax):
        for j, col in enumerate(row):
            col.set_xlabel('x')
            col.set_ylabel('u')

    lines, labels = fig.axes[-1].get_legend_handles_labels()
    fig.legend(lines, labels,
               bbox_to_anchor=(1.00, 0.96), ncol=1, framealpha=1,
               fontsize=8)

    c_template = 'c = %.3f'
    c_text = fig.text(0.87, 0.96, '')
    c_text.set_text(c_template % (c))

    fig.suptitle(f'Calculated by the {scheme} scheme', x=0.40, y=0.97,
                 fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'./outputs/u_{scheme}_{c:.3f}.png')
    plt.show()
    plt.close()


def animation(x, t, u, dt, c, scheme, ue):
    fig, ax = plt.subplots()
    plotu, = ax.plot([], [], color='cornflowerblue', linestyle='-', marker='o', label=f'u_{scheme}', animated=False)
    plotue, = ax.plot([], [], color='lightcoral', linestyle='-.', label='exact', animated=False)
    time_template = 'time = %.2fs'
    time_text = ax.text(0.88, 1.02, '', transform=ax.transAxes)
    c_template = 'c = %.3f'
    c_text = ax.text(0.88, 1.06, '', transform=ax.transAxes)
    c_text.set_text(c_template % (c))

    def init():
        ax.set_xlim(0, 3)
        ax.set_ylim(-1, 1)
        ax.set_xlabel('x', fontsize=12, fontweight='bold')
        ax.set_ylabel('u', fontsize=12, fontweight='bold')
        ax.grid(linestyle='-.')
        ax.set_title(f'Calculated by the {scheme} scheme', x=0.45, y=1.0,
                     fontsize=12, fontweight='bold')
        time_text.set_text('')
        ax.legend(loc='lower left')
        return plotu, plotue, time_text,

    def update(i):
        newx = x
        newu = u[i]
        newue = ue[i]
        plotu.set_data(newx, newu)
        plotue.set_data(newx, newue)
        time_text.set_text(time_template % (dt * i))
        return plotu, plotue, time_text,

    ani = FuncAnimation(fig, update, len(t) - 1,
                        init_func=init)
    ani.save(f'./outputs/{scheme}_{c:.3f}.gif', writer='ffmpeg', fps=60)
    plt.show()
    plt.close()


# 生成误差图
def plot_l2error(dxt, l2err, scheme, p, check, Theory):
    delta = l2err[0] / dxt[0] ** Theory
    ref = [d ** Theory * delta for d in dxt]
    plt.loglog(dxt, l2err, linewidth=3.0, linestyle='-', marker='s', label=scheme + f', p={p:.2f}')
    plt.xticks(np.array(dxt))
    plt.loglog(dxt, ref, linewidth=3.0, linestyle='-.', marker='o', label=f'Theory, p ={Theory}')
    plt.xlabel(rf'd{check}', fontsize=15)
    plt.ylabel(r'$\|u_h-u\|$', fontsize=15)
    # plt.xticks(np.array(dxt))
    plt.title(f'Verify the accuracy order of the {scheme} in {check}', fontsize=12, fontweight='bold')
    plt.legend()
    plt.savefig(f'./outputs/l2error_{scheme}_{check}.png')
    plt.show()
    plt.close()
