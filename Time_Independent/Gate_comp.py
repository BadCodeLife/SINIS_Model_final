import Calculation_Methods_SS as ss
import setup_file as su
import scipy as sp
import matplotlib.pyplot as plt
import copy

file = 'SteadyStateOut.h5'
key = 'current'

input = ss.create_datapoint_dict_steady(
    gate=su.rounded_linspace(0, 0.5, 4),
    bias=su.rounded_linspace(-5, 5, 1001),
    n_spread=3,
    temp=1e-2,
    leak=1e-4,
    charge_energy=1,
    super=True
)

# ss.check_tested_settings(file,key)

if __name__ == '__main__':
    Existing_Data = ss.fetch_data(file, key)
    # print (Existing_Data.to_string())

    Gate = input[ss.Gate_name]
    Bias = input[ss.Bias_name]

    desired_points = []
    for gate in Gate:
        for bias in Bias:
            data_point = copy.deepcopy(input)
            data_point[ss.Gate_name] = gate
            data_point[ss.Bias_name] = bias
            desired_points.append(data_point)

    points_to_calc = ss.non_existing_points(desired_points, Existing_Data)

    print 'number of points to calculate', len(points_to_calc)

    ss.calculate_points(points_to_calc, file, key, Existing_Data)
    Existing_Data = ss.fetch_data(file, key)

    fig, ax = plt.subplots(1)
    colour = 0
    for gate in Gate:
        line_data_set = copy.deepcopy(input)
        line_data_set[ss.Gate_name] = gate
        line_label = '%s$\,$e' % gate
        ss.current_line_plot(line_data_set, Existing_Data, ax, colour, '-', line_label)
        colour+=1




    ax.legend(title='Gate Charge')

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

    plt.rcParams.update({'font.size': 20})
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_ylabel('Current [ef]', fontsize=20)  # Y label
    ax.set_xlabel('Bias Voltage [eV/$\Delta$]', fontsize=20)  # X label




    plt.suptitle(Title)
    plt.show()

