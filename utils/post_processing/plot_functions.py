import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, rc
from IPython.display import HTML


t_day = 3600 * 24.0
t_year = 365 * t_day
figsize = (9, 8)


def timeseries(ot, ot_vmax):

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=figsize)

    axes[0].plot(ot["t"] / t_year, ot["tau"])
    axes[0].set_ylabel("tau [Pa]")

    axes[1].plot(ot["t"] / t_year, ot["theta"])
    axes[1].set_ylabel("state [s]")

    axes[2].plot(ot_vmax["t"] / t_year, ot_vmax["v"])
    axes[2].set_ylabel("max v [m/s]")
    axes[2].set_xlabel("time [yr]")
    axes[2].set_yscale("log")

    plt.tight_layout()
    plt.show()


def slip_profile(ox, warm_up=0, orientation="horizontal"):
    x_unique = ox["x"].unique()
    sort_inds = np.argsort(x_unique)
    x_unique = x_unique[sort_inds]
    t_vals = np.sort(ox["t"].unique())

    if warm_up > t_vals.max():
        print("Warm-up time > simulation time!")
        return

    ind_warmup = np.where(t_vals >= warm_up)[0][0]

    Nx = len(x_unique)
    Nt = len(t_vals) - 1

    slice = np.s_[Nx * ind_warmup:Nx * Nt]
    data_shape = (Nt - ind_warmup, Nx)
    x = ox["x"][slice].values.reshape(data_shape)[:, sort_inds]
    slip = ox["slip"][slice].values.reshape(data_shape)[:, sort_inds]
    v = ox["v"][slice].values.reshape(data_shape)[:, sort_inds]

    slip -= slip[0]
    t = t_vals[ind_warmup:-1]
    t -= t[0]

    fig = plt.figure(figsize=figsize)

    if orientation == "horizontal":
        CS = plt.contourf(x, slip, np.log10(v), levels=200, cmap="magma")
        plt.xlabel("position [m]")
        plt.ylabel("slip [m]")
        CB = plt.colorbar(CS, orientation="horizontal")
        CB.ax.set_title("slip rate [m/s]")
    elif orientation == "vertical":
        x -= x.min()
        CS = plt.contourf(slip, x * 1e-3, np.log10(v), levels=200, cmap="magma")
        plt.ylabel("depth [km]")
        plt.xlabel("slip [m]")
        plt.gca().invert_yaxis()
        CB = plt.colorbar(CS, orientation="horizontal")
        CB.ax.set_title("slip rate [m/s]")
    else:
        print("Keyword 'orientation=%s' not recognised" % orientation)
        plt.close()
        return
    plt.tight_layout()
    plt.show()


def animation_slip(ox, warm_up=0, orientation="horizontal"):

    x_unique = ox["x"].unique()
    sort_inds = np.argsort(x_unique)
    x_unique = x_unique[sort_inds]
    t_vals = np.sort(ox["t"].unique())

    if warm_up > t_vals.max():
        print("Warm-up time > simulation time!")
        return

    ind_warmup = np.where(t_vals >= warm_up)[0][0]

    Nx = len(x_unique)
    Nt = len(t_vals) - 1

    slice = np.s_[Nx * ind_warmup:Nx * Nt]
    data_shape = (Nt - ind_warmup, Nx)
    slip = ox["slip"][slice].values.reshape(data_shape)[:, sort_inds]
    v = ox["v"][slice].values.reshape(data_shape)[:, sort_inds]

    slip -= slip[0]
    t = t_vals[ind_warmup:-1]
    t -= t[0]

    plt.ioff()

    if orientation == "horizontal":
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=figsize)

        ax[0].set_xlim((x_unique.min(), x_unique.max()))
        ax[0].set_ylim((slip.min(), slip.max()))
        ax[0].set_ylabel("slip [m]")

        ax[1].set_xlim((x_unique.min(), x_unique.max()))
        ax[1].set_ylim((v.min(), v.max()))
        ax[1].set_yscale("log")
        ax[1].set_ylabel("slip rate [m/s]")
        ax[1].set_xlabel("position [m]")

        line1, = ax[0].plot([], [], lw=2)
        line2, = ax[1].plot([], [], lw=2)
        lines = (line1, line2)

        def init_plot():
            for line in lines:
                line.set_data([], [])
            return lines,

        def animate(i):
            line1.set_data(x_unique, slip[i])
            line2.set_data(x_unique, v[i])

    elif orientation == "vertical":
        x_unique -= x_unique.min()
        x_unique *= 1e-3
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=figsize)

        ax[0].set_ylim((x_unique.min(), x_unique.max()))
        ax[0].set_xlim((slip.min(), slip.max()))
        ax[0].set_ylabel("depth [km]")
        ax[0].set_xlabel("slip [m]")
        ax[0].invert_yaxis()
        ax[0].xaxis.tick_top()
        ax[0].xaxis.set_label_position("top")

        ax[1].set_ylim((x_unique.min(), x_unique.max()))
        ax[1].set_xlim((v.min(), v.max()))
        ax[1].set_xscale("log")
        ax[1].set_xlabel("slip rate [m/s]")
        ax[1].invert_yaxis()
        ax[1].xaxis.tick_top()
        ax[1].xaxis.set_label_position("top")

        line1, = ax[0].plot([], [], lw=2)
        line2, = ax[1].plot([], [], lw=2)
        lines = (line1, line2)

        def init_plot():
            for line in lines:
                line.set_data([], [])
            return lines,

        def animate(i):
            line1.set_data(slip[i], x_unique)
            line2.set_data(v[i], x_unique)

    anim = animation.FuncAnimation(fig, animate, init_func=init_plot,
                                   frames=len(t), interval=30, blit=True)
    plt.tight_layout()
    HTML(anim.to_html5_video())
    rc("animation", html="html5")
    plt.close()
    return anim

def plot_2D_slip_rate(ox, t_beg=None,t_end=None,V_pl=1e-9,Lvw=5,Lc=None,normalize=True,
                              cmap='RdBu_r',vmin=None,vmax=None,is_log=True):

    x_unique = ox["x"].unique()
    sort_inds = np.argsort(x_unique)
    x_unique = x_unique[sort_inds]
    t_vals = np.sort(ox["t"].unique())

    cdt = (t_vals >= t_beg) & (t_vals <= t_end)
    ind_beg, ind_end = np.where(cdt)[0][0], np.where(cdt)[0][-1]
    print (t_vals[ind_beg], t_vals[ind_end])

    Nx = len(x_unique)
    Nt = ind_end- ind_beg
    slice = np.s_[Nx* ind_beg:Nx * ind_end]
    data_shape = (Nt, Nx)
    print ('data_shape (Nt, Nx):', data_shape)
    x = ox["x"][slice].values.reshape(data_shape)[:, sort_inds]
    v = ox["v"][slice].values.reshape(data_shape)[:, sort_inds]
    t = ox["t"][slice].values.reshape(data_shape)[:, sort_inds]
    print (x.shape, v.shape, t.shape)
    x /= Lc
    t -= t[0]
    t /= t_year
    if normalize: v /= V_pl

    data = v
    if is_log: data = np.log10(v)
    if vmin==None: vmin= data.min()
    if vmax==None: vmax= data.max()
    print ('vmin, vmax: ', vmin, vmax)

    fig = plt.figure(figsize=(8, 4))
    ext = [np.amin(t), np.amax(t), np.amin(x), np.amax(x)]
    im = plt.imshow(data.T, extent=ext, aspect='auto',cmap=plt.cm.get_cmap(cmap,31),vmin=vmin,vmax=vmax)
    plt.xlabel("Time since warm up (yr)")
    plt.ylabel("Along fault (Lnuc)")
    plt.axhline(Lvw, c='snow',lw=0.5, ls=':')
    plt.axhline(-Lvw, c='snow',lw=0.5, ls=':')
    plt.ylim(-Lvw*1.5, Lvw*1.5)
    if normalize: plt.colorbar(im, shrink=0.7, label='Slip rate, log10(V/V_pl)')
    if not normalize: plt.colorbar(im, shrink=0.7, label='Slip rate, V (m/s)')
    plt.title("Warm up (yr) = " + str(t_beg/t_year))
###