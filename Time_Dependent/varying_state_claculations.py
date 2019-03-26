import setup_file as su
import scipy as sp
import pandas as pd
import multiprocessing as mp
import copy

Gate_Amp_name = 'Gate_Amplitude_(e)'
Gate_Function_name = 'Gate_Function'
Bias_Func_name = 'Bias_Function'  # Assumption at this point is that the variance is between 0 and 1
States_name = 'States Set'  # States of the island
Thermal_E_name = 'Thermal_Energy_(Delta)'
Leak_name = 'Leakage_(Delta)'
Super_Con_name = 'Superconducting'
Period_time_name = 'Time Period'
Period_number_name = 'Number of Periods'
Total_Step_name = 'Number of Steps'
Current_name = 'Current'


def create_data_entry(gate_amp,gate_func,bias_func,states,temp,leak,super,period,nu_periods,time_steps):
    dict =