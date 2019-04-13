import setup_file as su
import scipy as sp
import pandas as pd
import multiprocessing as mp
import copy

Gate_name = 'Gate_(e)'
Bias_name = 'Bias_(eV)'
States_name = 'States Set'
Thermal_E_name = 'Thermal_Energy_(Delta)'
Leak_name = 'Leakage_(Delta)'
Charge_E_name = 'Charging_Energy_(Delta)'  # charging energy of the island
Super_Con_name = 'Superconducting'
Current_name = 'Current'

# Method that makes a more formal definition of what the expected dictionary for a data point should be.
# This is an incomplete DataPoint as it does not include the expectation value of current 'Current'
def create_datapoint_dict_steady(gate,bias,n_set,temp,leak,charge_energy,super):
    dict = {
        Gate_name:gate,
        Bias_name:bias,
        States_name:n_set,
        Thermal_E_name:temp,
        Leak_name:leak,
        Charge_E_name: charge_energy,
        Super_Con_name:super
    }
    return dict


########################################################################################################################
# calculations

# calculates the
def create_matrix(gate,bias,n_set):
    dim = len(n_set)
    matrix = sp.zeros((dim,dim))
    for i in xrange(dim):
        if i<dim-1:
            matrix[i][i] = su.tunnel_on(gate,bias,n_set[i])
            matrix[i][i+1] = -su.tunnel_off(gate,bias,n_set[i+1])
        matrix[dim-1][i] = 1
    return matrix


def probability_calc(gate,bias,n_set):
    dim = len(n_set)
    y = sp.zeros(dim)
    y[dim-1]=1
    matrix = create_matrix(gate,bias,n_set)
    out = sp.linalg.solve(matrix,y)
    return out


def tunnel_vector(gate,bias,n_set):
    out = sp.array([])
    for n in n_set:
        rate = su.tunnel_flow_j1(gate_charge=gate,bias_potential=bias,extra_electrons=n)
        out = sp.append(out,rate)
    return out


def steady_state(gate,bias,n_set):
    p = probability_calc(gate,bias,n_set)
    t = tunnel_vector(gate,bias,n_set)
    return sp.dot(p,t)

########################################################################################################################

def fetch_data(data_file,file_key):
    try:
        fetched_data = pd.DataFrame(pd.read_hdf(data_file,file_key,mode='r'))
    except:
        fetched_data = pd.DataFrame(columns={Gate_name, Bias_name, States_name, Thermal_E_name,
                                             Leak_name, Super_Con_name, Current_name})
    return fetched_data

def data_query(data_point, existing_data):
    point_locations = existing_data[
        (existing_data[Gate_name] == data_point[Gate_name]) &
        (existing_data[Bias_name] == data_point[Bias_name]) &
        (existing_data[States_name] == data_point[States_name]) &
        (existing_data[Thermal_E_name] == data_point[Thermal_E_name]) &
        (existing_data[Leak_name] == data_point[Leak_name]) &
        (existing_data[Charge_E_name] == data_point[Charge_E_name]) &
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


def calculate_points(points_to_calc, data_file, file_key, existing_data):
    if len(points_to_calc) != 0:
        with mp.Manager() as manager:
            point = 0
            while point < len(points_to_calc):
                calculated_points = manager.list()
                processes = []
                for i in xrange(mp.cpu_count() - 2):
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


def current_point(point,calculated_points):
    su.thermal_energy = point[Thermal_E_name]
    su.leakage = point[Leak_name]
    su.super_conductor = point[Super_Con_name]
    su.unit_charge_energy = point[Charge_E_name]

    current = steady_state(gate=point[Gate_name],bias=point[Bias_name],n_set=point[States_name])

    complete_point = copy.deepcopy(point)
    complete_point[Current_name]=current

    calculated_points.append(complete_point)


colour_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                   '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5']


# This method will print a line of the plot. the line is going to be with the gate amplitude as the x-axis
def current_line_plot(line_data_set, existing_data, axis, colour, line_label):
    Bias = line_data_set[Bias_name]
    del line_data_set[Bias_name]

    Current_list = []
    for bias in Bias_name:
        line_data_set[Bias_name] = bias

        data_point_index = data_query(line_data_set, existing_data)[0]

        current = existing_data.at[data_point_index, 'Current']
        Current_list.append(current)

    axis.plot(Bias, Current_list, '.-', label=line_label, color=colour_sequence[colour])


# This method will print a table with all the settings that have been tested and are inside the datafile
def check_tested_settings(file,key):
    fetched_data = pd.DataFrame(pd.read_hdf(file, key, mode='r'))
    del fetched_data[Gate_name]
    del fetched_data[Current_name]
    wanted_info = [Bias_name, States_name, Thermal_E_name, Leak_name, Charge_E_name, Super_Con_name ]

    fetched_data = fetched_data.drop_duplicates(subset=wanted_info, keep='first')

    print(fetched_data.to_string())