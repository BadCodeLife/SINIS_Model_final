import setup_file as su
import scipy as sp
import matplotlib.pyplot as plt

def create_datapoint_dict_steady(gate,bias,n_set,temp,leak,super):
    dict = {
        'Gate':gate,
        'Bias':bias,
        'States Set':n_set,
        'Thermal Energy':temp,
        'Leakage':leak,
        'Superconducting':super
    }
    return dict
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


print su.dos.__name__
su.super_conductor = True
su.check_sup()
print su.dos.__name__
# gate = sp.linspace(0,1,15)
# n_set = sp.arange(-3,4)
# bias = 0
# current = []
# for g in gate:
#     current.append(steady_state(g,bias,n_set))
#
# plt.plot(gate,current)
# plt.show()