import matplotlib.pyplot as plt
import scipy as sp
import ss_calculation as ss


Bias = sp.linspace(-5,5,1001)
gate = 0.
n_set = sp.arange(-7,8)

current = []
i = 0
for bias in Bias:
    i+=1
    print i
    current.append(ss.steady_state(gate,bias,n_set))

plt.plot(Bias,current)
plt.show()