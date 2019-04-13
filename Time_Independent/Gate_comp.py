import Calculation_Methods_SS as ss
import setup_file as su
import scipy as sp
import matplotlib.pyplot as plt
import copy

file = 'SteadyStateOut.h5'
key = 'current'

input = ss.create_datapoint_dict_steady(
    gate=su.rounded_linspace(0, 1, 7),
    bias=su.rounded_linspace(-5, 5, 1001),
    n_set=sp.arange(-7, 8),
    temp=1e-2,
    leak=1e-4,
    charge_energy=1,
    super=True
)


Existing_Data = ss.fetch_data(file, key)

Gate = input[ss.Gate_name]
Bias = input[ss.Bias_name]

desired_points = []
for gate in Gate:
    for bias in Bias:
        data_point = copy.deepcopy(input)
        data_point[ss.Gate_name] = gate
        data_point[ss.Bias_name] = bias
        desired_points.append(data_point)

