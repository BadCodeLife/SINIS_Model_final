import scipy as sp
import scipy.linalg as la
import setup_file as su
import multiprocessing as mp
import pandas as pd
import copy

# print hasattr(b,"__len__")

Gate_Amp_name = 'Gate_Amplitude_(e)'  # Amplitude of gate charge
Gate_Func_name = 'Gate_Function'  # Function of gate charge change
Gate_Occ_Cent_name = 'Gate_Oscillation_Center_(e)'
Bias_Func_name = 'Bias_Function'  # function of how bias varies with time
Time_Period_name = 'Period_of_Oscillation_(RC)'  # period of oscillation for a single cycle of the transistor
Number_of_Periods_name = 'Number_of_Periods'  # number of periods of oscillation that are being calculated
Number_of_Steps_name = 'Number_of_Time_Steps'  # number of steps all the periods will be split into
States_name = 'States Range'
Thermal_E_name = 'Thermal_Energy_(Delta)'
Leak_name = 'Leakage_(Delta)'
Charge_E_name = 'Charging_Energy_(Delta)'  # charging energy of the island
Super_Con_name = 'Superconducting'
Current_name = 'Current_(ef)'


# Method that makes a more formal definition of what the expected dictionary for a data point should be.
# This is an incomplete DataPoint as it does not include the expectation value of current 'Current'
def create_data_point_dict(gate_amp, gate_func, gate_occ_cent, bias_function, period, number_of_periods,
                           time_steps, n_set, temp, leak, charge_energy, superconductor):
    dictionary = {
        Gate_Amp_name: gate_amp,
        Gate_Func_name: gate_func,
        Gate_Occ_Cent_name: gate_occ_cent,
        Bias_Func_name: bias_function,
        Time_Period_name: period,
        Number_of_Periods_name: number_of_periods,
        Number_of_Steps_name: time_steps,
        States_name: n_set,
        Thermal_E_name: temp,
        Leak_name: leak,
        Charge_E_name: charge_energy,
        Super_Con_name: superconductor
    }
    return dictionary


# Method that creates the array of times used for this calculation.
# This is under the assumption that for N steps, you want to use t_0 -> t_N-1
def time_array_calc(charging_energy, time_period, number_of_periods, number_of_steps):
    end_time = time_period * number_of_periods / charging_energy
    time_array = su.rounded_linspace(0, end_time, number_of_steps+1)
    time_array = time_array[:number_of_steps]
    return time_array


########################################################################################################################
# makes the matrix that is defined by the relation dP(n)/dt = gamma(+)P(n-1)+gamma(-)P(n+1)-gamma(o)P(n).
# the matrix is defiend by the normalised gate voltage and normalised bias voltage.
# n_set defines the values of n considered on the diagonal.
# matrix assumes that P(n+1) = 0 on the upper end of n_set and P(n-1) on the lower end of n_set.
# This is an old method that was used as a concept. For an improved version, see make_diff_matrix
def make_matrix(gate, bias, n_set):
    dim = len(n_set)
    matrix = sp.zeros((dim, dim))
    for i in xrange(dim):
        matrix[i, i] = -su.tunnel_total(gate, bias, n_set[i])
        if i < dim - 1:
            matrix[i, i + 1] = su.tunnel_off(gate, bias, n_set[i + 1])
        if i > 0:
            matrix[i, i - 1] = su.tunnel_on(gate, bias, n_set[i - 1])
    return matrix


# creates a vector that defines the overall tunnelling rate for junction 1 for all n states in n_set.
# This is an old version used as proof of concept. see make_tunnel_vector_v2
def make_tunnel_vector(gate, bias, n_set):
    dim = len(n_set)
    out = sp.array([])
    for i in xrange(dim):
        out = sp.append(out, su.tunnel_flow_j1(gate, bias, n_set[i]))
    return out


# creates a matrix of the reciprocal relationship in make_matrix, as well as appending bellow the vector from
# make_tunnel_vector. appended to the side it makes the vector (0,1) which allows the calculation to produce a matrix
# for which an exponent exists.
def expand_matrix(gate, bias, n_set):
    dim = len(n_set)
    matrix = make_matrix(gate, bias, n_set)
    vector = make_tunnel_vector(gate, bias, n_set).reshape((1, dim))
    zeros = sp.zeros((dim + 1, 1))

    A = sp.concatenate((matrix, vector), axis=0)
    A = sp.concatenate((A, zeros), axis=1)
    return A


# exponentiates the expanded matrix, creating a time evolution operator for the vector P.
# n_set defines the states considered during the operation cycle and time_array defines the time steps used to
# calculate it.
def exponential_matrix(gate_amp, gate_function, bias_func, n_set, time_array):
    time_step = time_array[1]
    U = None
    for time in time_array:
        gate = gate_function(gate_amp, time)
        bias = sp.float64(bias_func(time))

        A = time_step * make_expanded_diff_matrix(gate, bias, n_set)
        segment = la.expm(A)
        if U is None:
            U = segment
        else:
            U = sp.matmul(segment, U)
    return U


# finds the probability distribution at the start of the operation cycle.
# operates under the assumption that distribution at the start of the cycle is the same as at the end.
def probability_calculation(U, n_set):
    dim = len(n_set)
    u = U[:dim, :dim]

    i = sp.identity(dim)
    m = sp.subtract(u, i)
    w = la.svd(m)[2]
    v = w[len(w) - 1]
    vector = v / sp.sum(v)

    return vector

# creates an expanded probability vector from the one calculated at the start of the cycle. The final value in the
# array represents the current transferred so far in electrons.
def expand_probability(p):
    return sp.append(p, 0)


# finds the current over the number of cycles defined by the time_array.
def current_calc(gate_amp, gate_func, bias_func, n_set, time_array, number_of_periods):
    dim = len(n_set)
    U = exponential_matrix(gate_amp, gate_func, bias_func, n_set, time_array)
    probability = probability_calculation(U, n_set)
    p_tilde = expand_probability(probability)
    tunneling_vec = U[dim]
    current = sp.dot(tunneling_vec, p_tilde)

    return current / number_of_periods


# will find all tunnel events possible for a n_g, V combination for the states in n_set
# output array is two dimensional, with [n][1+,1-,2+,2-] as the locations
# used to cut down on calculation times.
def tunnel_array(gate, bias, n_set):
    output_array = sp.array([])
    for n in n_set:
        tunnel_events = sp.array([[su.tunnel_rate_1p(gate, bias, n),
                                   su.tunnel_rate_1n(gate, bias, n),
                                   su.tunnel_rate_2p(gate, bias, n),
                                   su.tunnel_rate_2n(gate, bias, n)]])
        if not output_array.any():
            output_array = tunnel_events
        else:
            output_array = sp.concatenate((output_array, tunnel_events), axis=0)
    return output_array


# makes the matrix that is defined by the relation dP(n)/dt = gamma(+)P(n-1)+gamma(-)P(n+1)-gamma(o)P(n).
# the matrix is defiend by the normalised gate voltage and normalised bias voltage.
# n_set defines the values of n considered on the diagonal.
# matrix assumes that P(n+1) = 0 on the upper end of n_set and P(n-1) on the lower end of n_set.
# This method uses a tunnelling array (t_array) to save on repeated calculations. This method saved on calculation times
# by roughly 60%.
def make_diff_matrix(t_array):
    dim = len(t_array)
    matrix = sp.zeros((dim, dim))
    for i in xrange(dim):
        matrix[i, i] = -(t_array[i, 0] + t_array[i, 1] + t_array[i, 2] + t_array[i, 3])
        if i < dim - 1:
            matrix[i, i + 1] = t_array[i + 1, 1] + t_array[i + 1, 3]
        if i > 0:
            matrix[i, i - 1] = t_array[i - 1, 0] + t_array[i - 1, 2]
    return matrix


# creates a vector that defines the overall tunnelling rate for junction 1 for all n states in n_set.
# This method saves on repeated calculations by using a t_array, containing all the relevent tunneling values.
def make_tunnel_vec_v2(t_array):
    dim = len(t_array)
    out = sp.array([])
    for i in xrange(dim):
        tunnel_current = t_array[i, 0] - t_array[i, 1]
        out = sp.append(out,tunnel_current)
    return out


# creates a matrix of the reciprocal relationship in make_matrix, as well as appending bellow the vector from
# make_tunnel_vector. appended to the side it makes the vector (0,1) which allows the calculation to produce a matrix
# for which an exponent exists.
# This method uses t_array to cut down on calculation time in make_diff_matrix and make_tunnel_vec_v2. This is
# because, most of the values are repeated between different values inside the recurrence relationship outlined in the
# thesis.
def make_expanded_diff_matrix(gate, bias, n_set):
    dim = len(n_set)
    t_array = tunnel_array(gate, bias, n_set)

    diff_mat = make_diff_matrix(t_array)
    tunnel_vec = make_tunnel_vec_v2(t_array).reshape((1, dim))
    zeros = sp.zeros((dim + 1, 1))

    A = sp.concatenate((diff_mat, tunnel_vec), axis=0)
    A = sp.concatenate((A, zeros), axis=1)

    return A

########################################################################################################################

# gathers existing data, as well as defining data if the files are empty.
def fetch_data(data_file, file_key):
    try:
        fetched_data = pd.DataFrame(pd.read_hdf(data_file, file_key, mode='r'))
    except:
        fetched_data = pd.DataFrame(columns={Gate_Amp_name, Gate_Func_name, Gate_Occ_Cent_name, Bias_Func_name,
                                             Time_Period_name, Number_of_Periods_name, Number_of_Steps_name,
                                             States_name, Thermal_E_name, Leak_name,Charge_E_name, Super_Con_name,
                                             Current_name})
    return fetched_data

# finds a data point within the existing data. This uses dictionary methods to confirm that the names are being
# used consistently. the method returns the index location of the data point within the existing data. This is in the
# form of the array.
# this is because the editor I used could help me keep track of the names since a lot of terms were being compared.
# this method is slightly slow when it comes to retrieving data, and this might be something worth
# improving in later iterations.
def data_query(data_point, existing_data):
    point_locations = existing_data[
        (existing_data[Gate_Amp_name] == data_point[Gate_Amp_name]) &
        (existing_data[Gate_Func_name] == data_point[Gate_Func_name].__name__) &
        (existing_data[Gate_Occ_Cent_name] == data_point[Gate_Occ_Cent_name]) &
        (existing_data[Bias_Func_name] == data_point[Bias_Func_name].__name__) &
        (existing_data[Time_Period_name] == data_point[Time_Period_name]) &
        (existing_data[Number_of_Periods_name] == data_point[Number_of_Periods_name]) &
        (existing_data[Number_of_Steps_name] == data_point[Number_of_Steps_name]) &
        (existing_data[States_name] == data_point[States_name]) &
        (existing_data[Thermal_E_name] == data_point[Thermal_E_name]) &
        (existing_data[Leak_name] == data_point[Leak_name]) &
        (existing_data[Charge_E_name] == data_point[Charge_E_name]) &
        (existing_data[Super_Con_name] == data_point[Super_Con_name])
        ].index
    return point_locations


# this methods confirms if one of the required data points is within the existing data.
# if multiple indices are returned it will notify that something has gone wrong leading to a data point being
# calculated twice.
def non_existing_points(data_points, existing_data):
    points_to_calc = list(data_points)
    for points in data_points:
        number_of_repeat_points = data_query(points, existing_data)
        if len(number_of_repeat_points) == 0:
            pass
        elif len(number_of_repeat_points) == 1:
            points_to_calc.remove(points)
        else:
            print('point has mulitple values, please check:%s')
            for x in points:
                print(x + ':' + str(points[x]))
            quit()
    return points_to_calc

# divides up the calculation of data points among cpu cores available. This leaves 2 cores unused. This is so the device
# can still be used while the calculations are being preformed. If this isn't desired, number_of_cores can be used to
# change the number of cores used during calculation.
# once all the cores have done the one data point assigned to them, they will save the existing data file before moving
# onto the next set of data points.
# information of how many data points are left will be provided.
# This is most likely not the most efficient method of preforming multiprocessing in this code.
def calculate_points(points_to_calc, data_file, file_key, existing_data, number_of_cores=(mp.cpu_count() - 2)):
    if len(points_to_calc) != 0:
        with mp.Manager() as manager:
            point = 0
            while point < len(points_to_calc):
                calculated_points = manager.list()
                processes = []
                for i in xrange(number_of_cores):
                    try:
                        p = mp.Process(target=current_point, args=(points_to_calc[point], calculated_points))
                        p.start()
                        processes.append(p)
                        point += 1
                        print point, '/', len(points_to_calc)
                    except:
                        pass

                for process in processes:
                    process.join()

                for points in calculated_points:
                    existing_data = existing_data.append(points, ignore_index=True)

                existing_data.to_hdf(data_file, key=file_key, mode='w')
    else:
        print('points exist for this line')

# this is the method that calls upon the calculation to find the current depending on the settings required by point.
# it appends the results onto an array held by calculate_points.
def current_point(point, calculated_points):
    su.super_conductor = point[Super_Con_name]

    su.gate_occ_center = point[Gate_Occ_Cent_name]
    su.thermal_energy = point[Thermal_E_name]
    su.leakage = point[Leak_name]

    su.unit_charge_energy = point[Charge_E_name]
    su.time_period = point[Time_Period_name] / point[Charge_E_name]
    su.frequency = 1 / su.time_period


    n_set = sp.arange(-point[States_name]-1,point[States_name]+1)


    time_array = time_array_calc(charging_energy=point[Charge_E_name],
                                 time_period=point[Time_Period_name],
                                 number_of_periods=point[Number_of_Periods_name],
                                 number_of_steps=point[Number_of_Steps_name])

    current =  current_calc(gate_amp=point[Gate_Amp_name],
                            gate_func=point[Gate_Func_name],
                            bias_func=point[Bias_Func_name],
                            n_set=n_set,
                            time_array=time_array,
                            number_of_periods=point[Number_of_Periods_name])

    completed_point = copy.deepcopy(point)
    completed_point[Gate_Func_name]=completed_point[Gate_Func_name].__name__
    completed_point[Bias_Func_name]=completed_point[Bias_Func_name].__name__
    completed_point[Current_name] = current

    calculated_points.append(completed_point)


colour_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                   '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5']

# This method will print a line of the plot. the line is going to be with the gate amplitude as the x-axis
def current_line_plot(line_data_set, existing_data, axis, colour, line_type,line_label):
    Gate_amplitudes = line_data_set[Gate_Amp_name]
    del line_data_set[Gate_Amp_name]

    Current_list = []
    for gate_amp in Gate_amplitudes:
        data_point = copy.deepcopy(line_data_set)
        data_point[Gate_Amp_name] = gate_amp

        data_point_index = data_query(data_point, existing_data)[0]

        current = existing_data.at[data_point_index, 'Current_(ef)']
        Current_list.append(current)

    axis.plot(Gate_amplitudes, Current_list, line_type, label=line_label, color=colour_sequence[colour])


# This method will print a table with all the settings that have been tested and are inside the datafile
def check_tested_settings(file,key):
    fetched_data = pd.DataFrame(pd.read_hdf(file, key, mode='r'))
    del fetched_data[Gate_Amp_name]
    del fetched_data[Current_name]
    wanted_info = [Gate_Func_name,Gate_Occ_Cent_name, Bias_Func_name, Time_Period_name, Number_of_Periods_name,
                   Number_of_Steps_name, States_name, Thermal_E_name, Leak_name, Charge_E_name, Super_Con_name ]

    fetched_data = fetched_data.drop_duplicates(subset=wanted_info, keep='first')

    print(fetched_data.to_string())

if __name__ == '__main__':

    # test = []
    #
    # test_point = create_data_point_dict(
    #     gate_amp=1.,
    #     gate_func=su.gate_curve,
    #     gate_occ_cent=0.5,
    #     bias_function=su.bias_unitary,
    #     period=1e-5,
    #     number_of_periods=1,
    #     time_steps=1400,
    #     n_set=3,
    #     temp=1e-2,
    #     leak=1e-4,
    #     charge_energy=1.,
    #     superconductor=True
    # )
    #
    # current_point(test_point,test)
    # for item in test[0]:
    #     print '%s:\t%s'%(item, test[0][item])

    t = sp.linspace(0, 1, 21)
    n = sp.arange(-3, 2)
    dim = len(n)
    U = exponential_matrix(1, su.gate_curve, su.bias_unitary, n, t[:20])

    p = probability_calculation(U, n)
    p = expand_probability(p)
    tun = U[dim]
    c = sp.dot(p, tun)
    print c