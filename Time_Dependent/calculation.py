import scipy as sp
import setup_file as su
import scipy.linalg as la


def make_matrix(gate, bias, n_set):
    dim = len(n_set)
    matrix = sp.zeros((dim, dim))
    for i in xrange(dim):
        matrix[i, i] = -su.tunnel_total(gate, bias, n_set[i])
        if i < dim - 1:            matrix[i, i + 1] = su.tunnel_off(gate, bias, n_set[i + 1])
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
    U = exponentiate(gate_amp, gate_func, bias_func, n_set, time_array)
    probability = probability_calculation(U, n_set)
    p_tilde = expand_probability(probability)
    tunneling_vec = U[dim]
    current = sp.dot(tunneling_vec, p_tilde)

    return current / number_of_periods