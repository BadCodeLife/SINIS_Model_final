import Calculation_Methods_SS as ss
import setup_file as su
import scipy as sp
from mpl_toolkits.mplot3d import axes3d     # This is required, even if not directly used, to recognise 3D environments
import matplotlib.pyplot as plt
import copy

file = 'SteadyStateOut.h5'
key = 'current'

input = ss.create_datapoint_dict_steady(
    gate=su.rounded_linspace(-5, 5, 201),
    bias=su.rounded_linspace(-5, 5, 201),
    n_spread=3,
    temp=1e-2,
    leak=1e-4,
    charge_energy=1,
    super=False
)


def Get_Data(input, existing_data):
    Gate = input[ss.Gate_name]
    Bias = input[ss.Bias_name]

    X, Y = sp.meshgrid(Bias, Gate)

    Z = sp.array([[]])
    for gate in Gate:
        z = sp.array([])
        for bias in Bias:
            data_point = copy.deepcopy(input)
            data_point[ss.Gate_name] = gate
            data_point[ss.Bias_name] = bias
            try:
                data_point_index = ss.data_query(data_point, existing_data)[0]
                current = existing_data.at[data_point_index, ss.Current_name]
                z = sp.append(z, current)
            except:
                print ss.data_query(data_point, existing_data)
                print data_point
                quit()

        try:
            Z = sp.vstack((Z, z))
        except:
            Z = z

    return X, Y, Z


def Wire_Plot(Bias, Gate, Current):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.rcParams.update({'font.size': 20})
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.plot_wireframe(Bias, Gate, Current)
    ax.set_ylabel('Gate Charge ($n_g$)', fontsize=20, labelpad=20)  # Y label
    ax.set_xlabel('Bias ($\Delta$/e)', fontsize=20, labelpad=20)  # X label
    ax.set_zlabel('Current/($\Delta$/eR)', fontsize=20, labelpad=20) # Z label
    ax.tick_params(pad=8)
    plt.show()


def Contourf_Plot(Bias, Gate, Current):
    fig, ax = plt.subplots(1)
    ax.contourf(Bias, Gate, Current)
    plt.show()


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

    Bias, Gate, Current = Get_Data(input, Existing_Data)

    Wire_Plot(Bias, Gate, Current)
