import Calculation_Methods_SS as ss
import setup_file as su
import scipy as sp
import matplotlib.pyplot as plt

file = 'SteadyStateOut.h5'
key = 'current'

input = ss.create_datapoint_dict_steady(
    gate=su.rounded_linspace(-1, 1, 7),
    bias=su.rounded_linspace(-5, 5, 1001),
    n_set=sp.arange(-7, 8),
    temp=1e-2,
    leak=1e-4,
    super=True
)


