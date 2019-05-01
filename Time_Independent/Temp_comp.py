import Calculation_Methods_SS as ss
import setup_file as su
import scipy as sp
import matplotlib.pyplot as plt
import copy

file = 'SteadyStateOut.h5'
key = 'current'

input = ss.create_datapoint_dict_steady(
    gate=0,
    bias=su.rounded_linspace(-5, 5, 1001),
    n_spread=3,
    temp=[0.003, 0.01, 0.1, 0.3,0.5, 1.],
    leak=1e-4,
    charge_energy=1,
    super=False
)

# ss.check_tested_settings(file,key)

if __name__ == '__main__':
    Existing_Data = ss.fetch_data(file, key)
    # print (Existing_Data.to_string())

    Temp = input[ss.Thermal_E_name]
    Bias = input[ss.Bias_name]

    desired_points = []
    for temp in Temp:
        for bias in Bias:
            data_point = copy.deepcopy(input)
            data_point[ss.Thermal_E_name] = temp
            data_point[ss.Bias_name] = bias
            desired_points.append(data_point)

    points_to_calc = ss.non_existing_points(desired_points, Existing_Data)

    print 'number of points to calculate', len(points_to_calc)

    ss.calculate_points(points_to_calc, file, key, Existing_Data)
    Existing_Data = ss.fetch_data(file, key)

    fig, ax = plt.subplots(1)
    colour = 0
    for temp in Temp:
        line_data_set = copy.deepcopy(input)
        line_data_set[ss.Thermal_E_name] = temp
        line_label = '%s$\,\Delta$' % temp
        ss.current_line_plot(line_data_set, Existing_Data, ax, colour, '-', line_label)
        colour += 1

    type = None
    if input[ss.Super_Con_name]:
        type = 'SINIS'
    else:
        type = 'NININ'

    Title = 'Steady state for different Gate Charges %s transistor \n' % type \
            + 'Charging Energy = %s $\Delta$, ' % input[ss.Charge_E_name] \
            + 'Thermal Energy = %s $\Delta$, ' % input[ss.Thermal_E_name] \
            + 'Leakage = %s $\Delta$, \n' % input[ss.Leak_name] \
            + 'Spread of states calculated is %s' % input[ss.States_name]

    plt.suptitle(Title)
    plt.rcParams.update({'font.size': 20})
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_ylabel('Current [e]', fontsize=20)  # Y label
    ax.set_xlabel('Bias Voltage [eV/$\Delta$]', fontsize=20)  # X label
    ax.legend(title='Thermal Energy')

    plt.show()
