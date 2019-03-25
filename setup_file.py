import scipy as sp
import scipy.integrate as int
from scipy import signal

# gate charge is measured in charge as a multiple of e
# bias potential is measured in delta
energy_gap = 1.  # normalised by itself
total_capacitance = 0.5  # in units of e^2/energy_gap
unit_charge_energy = 1. / (2. * total_capacitance)  # in units of energy_gap
tunnel_resistance = 1.  # normalised by itself
thermal_energy = 1e-2  # in units of energy_gap => temp = energy_gap/boltzman_const.
leakage = 1e-4  # in units of energy_gap
super_conductor = False  # define whether we are using SINIS or NININ

int_limit_low  = sp.inf
int_limit_high = sp.inf

frequency = 1.  # in units of 1/RC --- condition is that f<< 1/RC
time_period = 1.  # in units of RC --- condition is that T>> RC
gate_occ_center = 0.  # in no. of electorns of charge
number_of_periods = 3
total_time = time_period*number_of_periods
steps_per_period = 200



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

dos = dos_NIN
def check_sup():
    if super_conductor:
        dos = dos_SIN

def potential_diff_1p(gate_charge, bias_potent, extra_electrons):
    return (+2. * unit_charge_energy * (extra_electrons + gate_charge + 0.5)) - (bias_potent / 2.)


def potential_diff_2p(gate_charge, bias_potent, extra_electrons):
    return (+2. * unit_charge_energy * (extra_electrons + gate_charge + 0.5)) + (bias_potent / 2.)


def potential_diff_1n(gate_charge, bias_potent, extra_electrons):
    return (-2. * unit_charge_energy * (extra_electrons + gate_charge - 0.5)) + (bias_potent / 2.)


def potential_diff_2n(gate_charge, bias_potent, extra_electrons):
    return (-2. * unit_charge_energy * (extra_electrons + gate_charge - 0.5)) - (bias_potent / 2.)


# rate electrons tunnel on to the island from junction 1
def tunnel_rate_1p(gate_charge, bias_potential, extra_electrons):
    if not super_conductor:
        density_of_states = dos_NIN
    if super_conductor:
        density_of_states = dos_SIN
    f = lambda e: density_of_states(e) * fermi(e) \
                  * (1 - fermi(e - potential_diff_1p(gate_charge, bias_potential, extra_electrons)))
    return int.quad(f, -int_limit_low, +int_limit_high)[0]


# rate electrons tunnel on to the island from junction 2
def tunnel_rate_2p(gate_charge, bias_potential, extra_electrons):
    if not super_conductor:
        density_of_states = dos_NIN
    if super_conductor:
        density_of_states = dos_SIN
    f = lambda e: density_of_states(e) * fermi(e) \
                  * (1 - fermi(e - potential_diff_2p(gate_charge, bias_potential, extra_electrons)))
    return int.quad(f, -int_limit_low, +int_limit_high)[0]


# rate electrons tunnel off of the island from junction 1
def tunnel_rate_1n(gate_charge, bias_potential, extra_electrons):
    if not super_conductor:
        density_of_states = dos_NIN
    if super_conductor:
        density_of_states = dos_SIN
    f = lambda e: density_of_states(e) * (1 - fermi(e)) \
                  * fermi(e + potential_diff_1n(gate_charge, bias_potential, extra_electrons))
    return int.quad(f, -int_limit_low, +int_limit_high)[0]


# rate electrons tunnel off of the island from junction 2
def tunnel_rate_2n(gate_charge, bias_potential, extra_electrons):
    if not super_conductor:
        density_of_states = dos_NIN
    if super_conductor:
        density_of_states = dos_SIN
    f = lambda e: density_of_states(e) * (1 - fermi(e)) \
                  * fermi(e + potential_diff_2n(gate_charge, bias_potential, extra_electrons))
    return int.quad(f, -int_limit_low, +int_limit_high)[0]


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

# different gate functions for the purpose of


def gate_square(gA, t):
    return gA*(-signal.square(2*sp.pi*(frequency*t+0.25)))+gate_occ_center


def gate_sine(gA,t):
    return gA*sp.sin(2*sp.pi*t/time_period)+gate_occ_center

def gate_triag(gA,t):
    if t<=time_period/2.:
        out = gA*(-1.+4.*t*frequency)
    elif t>time_period/2.:
        out = gA*(1-4.*(t-time_period/2.)*frequency)
    else:
        print('what did you do!')
        quit()
    return out

def gate_triag_2(gA,t):
    if t<=time_period/4.:
        out =  -gA*(-4.*t*frequency)
    elif t<= time_period*(3./4.):
        out =  -gA*(-1.+4.*(t-time_period/4.)*frequency)
    elif t>time_period*(3./4.):
        out = -gA*(1-4*(t-time_period*(3./4.))*frequency)
    else:
        print('What did you do!')
        quit()
    return out

def gate_triag_3(gA,t):
    sum = 0
    for i in xrange(100):
        n = 2*i +1
        sum += ((pow(-1.,i))/(pow(n,2)))*sp.sin(2*n*sp.pi*frequency*t)
    out = gA*(8/pow(sp.pi,2))*sum + gate_occ_center
    return out

def gate_box(gA,t):
    sum = 0
    for i in xrange(1000):
        n=2*i+1
        sum += sp.sin(2*n*sp.pi*frequency*t)/n
    out = gA*(4./sp.pi)*sum+gate_occ_center
    return out
#-----------------------------------------------------------------------------------------------------------------------

# bias and gate function
def gate_curve(gA,t):
    return gA * sp.sin(2 * sp.pi * t / time_period) + gate_occ_center

def bias_curve(t):
    return 1 - sp.square(sp.cos(2 * sp.pi * t/time_period))

def gate_triangle(gA, t):
    sum = 0
    for i in xrange(100):
        n = 2 * i + 1
        sum += ((pow(-1., i)) / (pow(n, 2))) * sp.sin(2 * n * sp.pi *  t/time_period)
    out = gA * (8 / pow(sp.pi, 2)) * sum
    return out

def gate_triangle_v2(gA, t):
    sum = 0
    for i in xrange(100):
        n = 2 * i + 1
        sum += ((pow(-1., i)) / (pow(n, 2))) * sp.sin(2 * n * sp.pi * t/time_period)
    out = gA * (8 / pow(sp.pi, 2)) * sum + gate_occ_center
    return out

def bias_triangle(t):
    sum = 0
    for i in xrange(100):
        n = 2 * i + 1
        sum += ((pow(-1., i)) / (pow(n, 2))) * sp.sin((2 * n * sp.pi * (2*t-0.75)/time_period))
    out = 0.5*(1-(8 / pow(sp.pi, 2)) * sum)
    return out

def bias_triangle_v2(t):
    sum = 0
    for i in xrange(100):
        n = 2 * i + 1
        sum += ((pow(-1., i)) / (pow(n, 2))) * sp.sin((2 * n * sp.pi * (2*t/time_period-0.25)))
    out = 0.5*(1-(8 / pow(sp.pi, 2)) * sum)
    return out

def bias_triangle_v3(t):
    sum = 0
    for i in xrange(100):
        n = 2 * i + 1
        sum += ((pow(-1., i)) / (pow(n, 2))) * sp.sin((2 * n * sp.pi * (2*t/time_period-0.75)))
    out = 0.5*(1-(8 / pow(sp.pi, 2)) * sum)
    return out

def gate_cos(gA,t):
    return gA*sp.cos(2*sp.pi*(t/time_period))+gate_occ_center

def bias_cosh(t):
    scale = 10
    max = sp.cosh(scale)
    if(t<0):
        return 0
    elif(t<time_period/2):
        return sp.cosh(2*scale*(2*t/time_period-0.5))/max
    elif(t>time_period/2):
        while(t>time_period/2):
            t -= time_period/2
        return sp.cosh(2*scale*(2*t/time_period-0.5))/max


def bias_box_v2(t):
    sum = 0
    for i in xrange(1000):
        n=2*i+1
        sum += sp.sin(2*n*sp.pi*(2*t/time_period+0.25))/n
    out = 0.5*(1+(4./sp.pi)*sum)
    return out


def bias_unitary(t):
    return 1.0