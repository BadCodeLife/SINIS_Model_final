import matplotlib.pyplot as plt
import scipy as sp
import Calculation_Methods_SS as ss
import setup_file as su


Bias = sp.linspace(-5,5,101)
gate = 0.
n_set = sp.arange(-7,8)

su.super_conductor = True

current = []
i = 0
for bias in Bias:
    i+=1
    print i
    current.append(ss.steady_state(gate,bias,n_set))

plt.plot(Bias,current)
plt.show()