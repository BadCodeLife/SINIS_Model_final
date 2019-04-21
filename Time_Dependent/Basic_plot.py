import Calculation_Methods_RF as cm
import setup_file as su
import scipy as sp
import matplotlib.pyplot as plt
import copy

file = 'TimeDepOut.h5'
key = 'current'
Log = False

input = cm.create_data_point_dict(
    gate_amp=su.rounded_linspace(0.1, 0.5, 41 ),
    gate_func=su.gate_curve,
    gate_occ_cent=-0.5,
    bias_function=su.bias_unitary,
    period=1e-3,
    number_of_periods=1,
    time_steps=150,
    n_set=4,
    temp=1e-5,
    leak=1e-10,
    charge_energy=1.,
    superconductor=True
)
if __name__ == '__main__':
    Existing_Data = cm.fetch_data(file, key)


    Gate_Amplitudes = input[cm.Gate_Amp_name]

    desired_points = []
    for gate_amp in Gate_Amplitudes:
        data_point = copy.deepcopy(input)
        data_point[cm.Gate_Amp_name] = gate_amp
        desired_points.append(data_point)

    points_to_calc = cm.non_existing_points(desired_points, Existing_Data)

    print 'number of points to calculate', len(points_to_calc)

    cm.calculate_points(points_to_calc, file, key, Existing_Data)
    Existing_Data = cm.fetch_data(file, key)

    fig, ax = plt.subplots(1)
    colour = 0
    line_data_set = copy.deepcopy(input)
    line_label = None
    cm.current_line_plot(line_data_set, Existing_Data, ax, colour, line_label)

    type = None
    if super:
        type = 'SINIS'
    else:
        type = 'NININ'

    n_high = input[cm.States_name]+1
    n_low = -input[cm.States_name]

    if Log:
        ax[0].semilogy()
        ax[0].set_ylim(0.9995, 1.0005)

    Title = 'Basic 2 Variable plot for a %s transistor \n' % type \
            + 'Gate Function = %s, ' % input[cm.Gate_Func_name].__name__ \
            + 'Gate Oscillation center = %s e, '%input[cm.Gate_Occ_Cent_name] \
            + 'Bias Function = %s, ' % input[cm.Bias_Func_name].__name__ \
            + 'Charging Energy = %s $\Delta$, ' % input[cm.Charge_E_name] \
            + 'Thermal Energy = %s $\Delta$, ' % input[cm.Thermal_E_name] \
            + 'Leakage = %s $\Delta$, \n' % input[cm.Leak_name] \
            + 'Time Period = %s RC, ' % input[cm.Time_Period_name] \
            + 'Number of Periods = %s, ' % input[cm.Number_of_Periods_name] \
            + 'Number of Time Steps = %s, ' % input[cm.Number_of_Steps_name] \
            + 'States considered are from n=%s to %s' % (n_low, n_high)

    plt.suptitle(Title)
    plt.rcParams.update({'font.size': 20})
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_ylabel('Current [ef]', fontsize=20)  # Y label
    ax.set_xlabel('Gate Charge Amplitude [e]', fontsize=20)  # X label
    plt.show()

    print Title