import scipy as sp
import setup_file as su
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm

def make_ticklabels_invisible(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)


def steady_state_log(ax):
    Bresolution = 801
    Gresolution = 801
    Gate = sp.linspace(-4, 4, Gresolution)

    File_name = 'CurrentSINIS_HighRes.txt'
    Current = sp.loadtxt(File_name)

    sp.asarray(Current)
    Current = sp.diff(Current)
    Current = Current.T
    Bias = sp.linspace(-8, 8, Bresolution - 1)
    ax.imshow(Current, cmap='jet', extent=[Gate.min(), Gate.max(),Bias.min(), Bias.max()], norm=LogNorm(), aspect='auto')



def plot_trajectories(ax1,ax2,ax3,gate,bias):
    su.time_period = 1
    su.frequency = 1
    T = sp.linspace(0,1,10001)
    su.gate_occ_center = 0.5
    gA = 0.5
    g_line = []
    b_line = []
    for t in T:
        g_line.append(gate(gA,t))
        b_line.append(bias(t))
    ax2.plot(T,g_line)
    ax3.plot(T,b_line)
    ax1.plot(g_line,b_line, color ='#ffbb78')




def map_trajectory(gate_function,bias_function):
    gs1 = GridSpec(3, 4)
    ax1 = plt.subplot(gs1[:-1, :])
    ax2 = plt.subplot(gs1[-1, :-2])
    ax3 = plt.subplot(gs1[-1, -2:])
    gs1.update(left=0.15, right=0.85, top=0.95, hspace=0.4, wspace=0.4)
    plt.rcParams.update({'font.size': 20})

    steady_state_log(ax1)
    ax1.set_xlabel('Gate Charge ($n_g$)', fontsize=20)  # Y label
    ax1.set_ylabel('Bias ($\Delta$/e)', fontsize=20)  # X label
    ax1.set_xlim(-0.5,1.5)
    ax1.set_ylim(-4,4)
    ax1.tick_params(axis='both', which='major', labelsize=20)


    plot_trajectories(ax1,ax2,ax3,gate_function,bias_function)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.set_xlabel('Time ($\\tau$)', fontsize=20)  # Y label
    ax2.set_ylabel('Gate Charge ($n_g$)', fontsize=20)

    ax3.tick_params(axis='both', which='major', labelsize=20)
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.tick_right()
    ax3.set_xlabel('Time ($\\tau$)', fontsize=20)  # Y label
    ax3.set_ylabel('Bias ($\Delta$/e)', fontsize=20)  # X label

# fig = plt.figure()
# map_trajectory(su.gate_curve,su.bias_unitary)
# plt.show()

def plot_stuff(ax,gate,bias,colour):
    su.time_period = 1
    su.frequency = 1
    T = sp.linspace(0,1,10001)
    su.gate_occ_center = 0.5
    gA = 0.5
    g_line = []
    b_line = []
    for t in T:
        g_line.append(gate(gA,t))
        b_line.append(bias(t))
    ax.plot(g_line,b_line,color=colour)

fig,ax = plt.subplots(1)
plt.rcParams.update({'font.size': 20})
colour_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                   '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5']
steady_state_log(ax)
plot_stuff(ax,su.gate_curve,su.bias_curve_000,colour_sequence[2])
plot_stuff(ax,su.gate_curve,su.bias_curve_050,colour_sequence[1])
plot_stuff(ax,su.gate_curve,su.bias_unitary,colour_sequence[0])

ax.set_xlabel('Gate Charge ($n_g$)', fontsize=20)  # Y label
ax.set_ylabel('Bias ($\Delta$/e)', fontsize=20)  # X label
ax.set_xlim(-0.5, 1.5)
ax.set_ylim(-4, 4)
ax.tick_params(axis='both', which='major', labelsize=20)
plt.show()
