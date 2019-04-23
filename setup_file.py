import scipy as sp
import scipy.integrate as integrate
from scipy import signal
from math import floor,log10

# gate charge is measured in charge as a multiple of e
# bias potential is measured in delta
energy_gap = 1.  # normalised by itself
unit_charge_energy = 1.  # in units of energy_gap and \propto 1/Csum   || due to normalisations: Ec =  1/Csum
tunnel_resistance = 1.  # normalised by itself
thermal_energy = 1e-2  # in units of energy_gap => temp = energy_gap/boltzman_const.
leakage = 1e-4  # in units of energy_gap
super_conductor = False  # define whether we are using SINIS or NININ

int_limit_low  = sp.inf
int_limit_high = sp.inf

frequency = 1.  # in units of 1/RC --- condition is that f<< 1/RC
time_period = 1.  # in units of RC --- condition is that T>> RC
gate_occ_center = 0.  # in no. of electrons of charge
number_of_periods = 3  # The code is designed to calculate the average current over n periods. 3 is the default
total_time = time_period*number_of_periods  # Technically this is not needed. RC is the unit measurement of time,
                                            # therefore current in ef only needs to be divided by number of periods.
number_of_steps = 1200  # Number of timesteps the total calculation is split into



def fermi(e):
    return 1 / (1 + sp.exp(e / thermal_energy))

def dos_NIN(e):
    return 1.

def dos_SIN(e):
    return sp.absolute(
        sp.real(
            (e + 1.j * leakage) / sp.sqrt(sp.square(e + 1.j * leakage) - 1.)
        )
    )



def potential_diff_1p(gate_charge, bias_potent, extra_electrons):
    return (+2. * unit_charge_energy * (extra_electrons + gate_charge + 0.5)) - (bias_potent / 2.)


def potential_diff_2p(gate_charge, bias_potent, extra_electrons):
    return (+2. * unit_charge_energy * (extra_electrons + gate_charge + 0.5)) + (bias_potent / 2.)


def potential_diff_1n(gate_charge, bias_potent, extra_electrons):
    return (-2. * unit_charge_energy * (extra_electrons + gate_charge - 0.5)) + (bias_potent / 2.)


def potential_diff_2n(gate_charge, bias_potent, extra_electrons):
    return (-2. * unit_charge_energy * (extra_electrons + gate_charge - 0.5)) - (bias_potent / 2.)


# rate electrons tunnel on to the island from junction 1
# use of density of states used depends on external boolean defining super_conductor
def tunnel_rate_1p(gate_charge, bias_potential, extra_electrons):
    if not super_conductor:
        density_of_states = dos_NIN
    if super_conductor:
        density_of_states = dos_SIN
    f = lambda e: density_of_states(e) * fermi(e) \
                  * (1 - fermi(e - potential_diff_1p(gate_charge, bias_potential, extra_electrons)))
    return integrate.quad(f, -int_limit_low, +int_limit_high)[0]


# rate electrons tunnel on to the island from junction 2
# use of density of states used depends on external boolean defining super_conductor
def tunnel_rate_2p(gate_charge, bias_potential, extra_electrons):
    if not super_conductor:
        density_of_states = dos_NIN
    if super_conductor:
        density_of_states = dos_SIN
    f = lambda e: density_of_states(e) * fermi(e) \
                  * (1 - fermi(e - potential_diff_2p(gate_charge, bias_potential, extra_electrons)))
    return integrate.quad(f, -int_limit_low, +int_limit_high)[0]


# rate electrons tunnel off of the island from junction 1
# use of density of states used depends on external boolean defining super_conductor
def tunnel_rate_1n(gate_charge, bias_potential, extra_electrons):
    if not super_conductor:
        density_of_states = dos_NIN
    if super_conductor:
        density_of_states = dos_SIN
    f = lambda e: density_of_states(e) * (1 - fermi(e)) \
                  * fermi(e + potential_diff_1n(gate_charge, bias_potential, extra_electrons))
    return integrate.quad(f, -int_limit_low, +int_limit_high)[0]


# rate electrons tunnel off of the island from junction 2
# use of density of states used depends on external boolean defining super_conductor
def tunnel_rate_2n(gate_charge, bias_potential, extra_electrons):
    if not super_conductor:
        density_of_states = dos_NIN
    if super_conductor:
        density_of_states = dos_SIN
    f = lambda e: density_of_states(e) * (1 - fermi(e)) \
                  * fermi(e + potential_diff_2n(gate_charge, bias_potential, extra_electrons))
    return integrate.quad(f, -int_limit_low, +int_limit_high)[0]


# total rate of electrons tunneling on to the island
def tunnel_on(gate_charge, bias_potential, extra_electrons):
    out = tunnel_rate_1p(gate_charge, bias_potential, extra_electrons) \
          + tunnel_rate_2p(gate_charge, bias_potential, extra_electrons)
    return out


# total rate of electrons tunneling off of the island
def tunnel_off(gate_charge, bias_potential, extra_electrons):
    out = tunnel_rate_1n(gate_charge, bias_potential, extra_electrons) \
          + tunnel_rate_2n(gate_charge, bias_potential, extra_electrons)
    return out


# total rate of tunneling across the island
def tunnel_total(gate_charge, bias_potential, extra_electrons):
    out = tunnel_rate_1p(gate_charge, bias_potential, extra_electrons)   \
          + tunnel_rate_2p(gate_charge, bias_potential, extra_electrons) \
          + tunnel_rate_1n(gate_charge, bias_potential, extra_electrons) \
          + tunnel_rate_2n(gate_charge, bias_potential, extra_electrons)
    return out


# total flow in the direction of the island across junction 1
def tunnel_flow_j1(gate_charge, bias_potential, extra_electrons):
    out = tunnel_rate_1p(gate_charge, bias_potential, extra_electrons) \
          - tunnel_rate_1n(gate_charge, bias_potential, extra_electrons)
    return out

#-----------------------------------------------------------------------------------------------------------------------
# Methods that allow for removing computational error from input parameters and in the case of time dependent methods,
# time steps. This is mainly because the gate value is on of the distinguishing

def round_sig(x, sig=7):
    try:
        return round(x, sig-int(floor(log10(abs(x))))-1)
    except:
        return 0


def round_array_sig(x,sig=7):
    out = sp.array([])
    for i in x:
        out = sp.append(out,round_sig(i,sig))
    return out


def rounded_linspace(start,end,steps,sig=7):
    array = sp.linspace(start,end,steps)
    out = round_array_sig(array,sig)
    return out


#-----------------------------------------------------------------------------------------------------------------------

# bias and gate function

def gate_curve(gA,t):
    return gA * sp.sin(2 * sp.pi * t / time_period) + gate_occ_center


def bias_curve_000(t):
    return 1 - sp.square(sp.cos(2 * sp.pi * t/time_period))


def bias_curve_050(t):
    return 1 - 0.5*sp.square(sp.cos(2 * sp.pi * t / time_period))


def bias_curve_080(t):
    return 1 - 0.2*sp.square(sp.cos(2 * sp.pi * t / time_period))


def bias_curve_090(t):
    return 1 - 0.1*sp.square(sp.cos(2 * sp.pi * t / time_period))


def gate_triangle(gA, t):
    sum = 0
    for i in xrange(100):
        n = 2 * i + 1
        sum += ((pow(-1., i)) / (pow(n, 2))) * sp.sin(2 * n * sp.pi * t/time_period)
    out = gA * (8 / pow(sp.pi, 2)) * sum + gate_occ_center
    return out


def bias_triangle_000(t):
    sum = 0
    for i in xrange(100):
        n = 2 * i + 1
        sum += ((pow(-1., i)) / (pow(n, 2))) * sp.sin((2 * n * sp.pi * (2*t/time_period-0.75)))
    out = 1-(0.5+0.5*(8 / pow(sp.pi, 2)) * sum)
    return out

def gate_box(gA,t):
    sum = 0
    for i in xrange(1000):
        n=2*i+1
        sum += sp.sin(n*(2*sp.pi*(t/time_period)))/n
    out = gA*(4./sp.pi)*sum + gate_occ_center
    return out

def bias_box_000(t):
    sum = 0
    for i in xrange(1000):
        n = 2 * i + 1
        sum += sp.sin(n * (2 * sp.pi * (2 * t / time_period) + sp.pi/4)) / n
    form = (4. / sp.pi) * sum
    out = 1-(0.5+0.5*form)
    return out


def bias_unitary(t):
    return 1.0