import scipy as sp
import scipy.linalg as la
import setup_file as su
import multiprocessing as mp
import pandas as pd
import copy

# print hasattr(b,"__len__")

Gate_Amp_name = 'Gate_Amplitude_(e)'  # Amplitude of gate charge
Gate_Func_name = 'Gate_Function'  # Function of gate charge change
Bias_Func_name = 'Bias_Function'  # function of how bias varies with time
Time_Period_name = 'Period_of_Oscillation_(RC)'  # period of oscillation for a single cycle of the transistor
Number_of_Periods_name = 'Number_of_Periods'  # number of periods of oscillation that are being calculated
Number_of_Steps_name = 'Number_of_Time_Steps'  # number of steps all the periods will be split into
States_name = 'States Set'
Thermal_E_name = 'Thermal_Energy_(Delta)'
Leak_name = 'Leakage_(Delta)'
Super_Con_name = 'Superconducting'
Current_name = 'Current'


# Method that makes a more formal definition of what the expected dictionary for a data point should be.
# This is an incomplete DataPoint as it does not include the expectation value of current 'Current'
def create_data_point_dict_steady(gate_amp, gate_func, bias_variance, bias_function, period, number_of_periods,
                                  time_steps, n_set, temp, leak, super):
    dictionary = {
        Gate_Amp_name: gate_amp,
        Gate_Func_name: gate_func,
        Bias_Func_name: bias_function,
        Time_Period_name: period,
        Number_of_Periods_name: number_of_periods,
        Number_of_Steps_name: time_steps,
        States_name: n_set,
        Thermal_E_name: temp,
        Leak_name: leak,
        Super_Con_name: super
    }
    return dictionary


# Method that creates the array of times used for this calculation.
# This is under the assumption that for N steps, you want to use t_0 -> t_N-1
def time_array_calc(time_period, number_of_periods, number_of_steps):
    end_time = time_period * number_of_periods
    time_array = su.rounded_linspace(0, end_time, number_of_steps)
    time_array = time_array[:number_of_steps + 1]
    return time_array


########################################################################################################################

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


def make_tunnel_vector(gate, bias, n_set):
    dim = len(n_set)
    out = sp.array([])
    for i in xrange(dim):
        out = sp.append(out, su.tunnel_flow_j1(gate, bias, n_set[i]))
    return out


def expand_matrix(gate, bias, n_set):
    dim = len(n_set)
    matrix = make_matrix(gate, bias, n_set)
    vector = make_tunnel_vector(gate, bias, n_set).reshape((1, dim))
    zeros = sp.zeros((dim + 1, 1))

    A = sp.concatenate((matrix, vector), axis=0)
    A = sp.concatenate((A, zeros), axis=1)
    return A


def exponentiate(gate_amp, gate_function, bias_func, n_set, time_array):
    time_step = time_array[1]
    U = None
    for time in time_array:
        gate = gate_function(gate_amp, time)
        bias = sp.float64(bias_func(time))

        A = time_step * expand_matrix(gate, bias, n_set)
        segment = la.expm(A)
        if U is None:
            U = segment
        else:
            U = sp.matmul(segment, U)
    return U


def probability_calculation(U, n_set):
    dim = len(n_set)
    u = U[:dim, :dim]

    i = sp.identity(dim)
    m = sp.subtract(u, i)
    w = la.svd(m)[2]
    v = w[len(w) - 1]
    vector = v / sp.sum(v)

    return vector


def expand_probability(p):
    return sp.append(p, 0)


def current_calc(gate_amp, gate_func, bias_func, n_set, time_array, number_of_periods):
    dim = len(n_set)
    U = exponentiate(gate_amp, gate_func,bias_func, n_set, time_array)
    probability = probability_calculation(U, n_set)
    p_tilde = expand_probability(probability)
    tunneling_vec = U[dim]
    current = sp.dot(tunneling_vec, p_tilde)

    return current / number_of_periods


########################################################################################################################

def fetch_data(data_file,file_key):
    try:
        fetched_data = pd.DataFrame(pd.read_hdf(data_file,file_key,mode='r'))
    except:
        fetched_data = pd.DataFrame(columns={Gate_Amp_name,Gate_Func_name, Bias_Func_name,
                                             Time_Period_name, Number_of_Periods_name, Number_of_Steps_name,
                                             States_name, Thermal_E_name, Leak_name, Super_Con_name, Current_name})
    return fetched_data


def data_query(data_point, existing_data):
    point_locations = existing_data[
        (existing_data[Gate_Amp_name] == data_point[Gate_Amp_name]) &
        (existing_data[Gate_Func_name] == data_point[Gate_Func_name]) &
        (existing_data[Bias_Func_name] == data_point[Bias_Func_name]) &
        (existing_data[Time_Period_name] == data_point[Time_Period_name]) &
        (existing_data[Number_of_Periods_name] == data_point[Number_of_Periods_name]) &
        (existing_data[Number_of_Steps_name] == data_point[Number_of_Steps_name]) &
        (existing_data[States_name] == data_point[States_name]) &
        (existing_data[Thermal_E_name] == data_point[Thermal_E_name]) &
        (existing_data[Leak_name] == data_point[Leak_name]) &
        (existing_data[Super_Con_name] == data_point[Super_Con_name])
     ]
    return point_locations


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


def calculate_points(points_to_calc, data_file, file_key, existing_data,number_of_cores=(mp.cpu_count()-2)):
    if len(points_to_calc) != 0:
        with mp.Manager() as manager:
            point = 0
            while point < len(points_to_calc):
                calculated_points = manager.list()
                processes = []
                for i in xrange(number_of_cores):
                    try:
                        p = mp.Process(target=current_point, args=(points_to_calc[point],calculated_points))
                        p.start()
                        processes.append(p)
                        point += 1
                        print point ,'/',len(points_to_calc)
                    except:
                        pass

                for process in processes:
                    process.join()

                for points in calculated_points:
                    existing_data = existing_data.append(points, ignore_index=True)

                existing_data.to_hdf(data_file, key=file_key, mode='w')
    else:
        print('points exist for this line')


def current_point(point, calculated_points):
    su.thermal_energy = point[Thermal_E_name]
    su.thermal_energy = point[Thermal_E_name]
    su.leakage = point[Leak_name]
    su.super_conductor = point[Super_Con_name]

    time_array = time_array_calc(time_period=point[Time_Period_name],
                                 number_of_periods=point[Number_of_Periods_name],
                                 number_of_steps=point[Number_of_Steps_name])

    current = current_calc(gate_amp=point[Gate_Amp_name],
                           gate_func=point[Gate_Func_name],
                           bias_func=point[Bias_Func_name],
                           n_set=point[States_name],
                           time_array=time_array,
                           number_of_periods=point[Number_of_Periods_name])


    completed_point = copy.deepcopy(point)
    completed_point[Current_name] = current

    calculated_points.append(completed_point)