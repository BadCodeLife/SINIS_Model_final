# a = {'a':1,'b':2}
# print a
# a['c']=3
# print a
import scipy as sp
import setup_file as su
import matplotlib.pyplot as plt
import pandas as pd


def test_1():
    time_period = 1

    def gate_triangle(gA, t):
        sum = 0
        for i in xrange(100):
            n = 2 * i + 1
            sum += ((pow(-1., i)) / (pow(n, 2))) * sp.sin(2 * n * sp.pi * t / time_period)
        out = gA * (8 / pow(sp.pi, 2)) * sum
        return out

    def bias_triangle_000(t):
        sum = 0
        for i in xrange(100):
            n = 2 * i + 1
            sum += ((pow(-1., i)) / (pow(n, 2))) * sp.sin((2 * n * sp.pi * (2 * t / time_period - 0.75)))
        out = 1 - (0.5 + 0.5 * (8 / pow(sp.pi, 2)) * sum)
        return out

    times = sp.linspace(0, 1, 101)
    plt.plot(gate_triangle(1, times), bias_triangle_000(times))
    plt.show()


def test_dataframes():
    dataframe = pd.DataFrame(columns={'a', 'b', 'c'})

    def input(n):
        return {'a': n, 'b': n+1, 'c': n+2}

    for n in xrange(5):
        dataframe = dataframe.append(input(n),ignore_index=True)

    print dataframe

    for index,row in dataframe.iterrows():
        print row['a'], row['b'], row['c']


def basic_test():
    print '%s %s'%(1,2)


def test_energy():
    e = sp.linspace(-2,2,1001)
    doe = su.dos_SIN
    plt.plot(e,doe(e))
    plt.show()


