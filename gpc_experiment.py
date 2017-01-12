#%%

import numpy
import scipy
import matplotlib.pyplot as plt

import os
import sys

sys.path.append('/home/bbales2/gpc')

import gpc

ts = numpy.linspace(0, 25.0, 200)

def f(s, t, a, b, d, g):
    x, y = s

    #y, yp = s
    dxdt = a * x - b * x * y
    dydt = d * x * y - g * y
    #dxdt = yp
    #dydt = -y
    return numpy.array([dxdt, dydt])

def func(x0, y0, a, b, d, g):
    y = numpy.log(scipy.integrate.odeint(f, [x0, y0], ts, (a, b, d, g)))

    return y

out = func(1.5, 2.0, 2.0 / 3.0, 4.0 / 3.0, 1.0, 1.0)

plt.plot(ts, numpy.exp(out[:, 0]), 'r')
plt.plot(ts, numpy.exp(out[:, 1]), 'b')
plt.legend(['rabbits', 'wolves'])
plt.show()
#%%

reload(gpc)

x0, y0, a, b, d, g = 1.5, 2.0, 2.0 / 3.0, 4.0 / 3.0, 1.0, 1.0

func2 = lambda x : func(x, y0, a, b, d, g)

hd = gpc.GPC(5, func2, [('n', (1.5, 0.25), 5)])

covar = hd.covar()

var = numpy.zeros((len(ts), 2))
for i in range(len(ts)):
    for j in range(2):
        var[i, j] = covar[i, j, i, j]
#%%
#plt.plot(ts, numpy.sqrt(var[:, 0]), 'r')
#plt.plot(ts, numpy.sqrt(var[:, 1]), 'b')
plt.plot(ts, numpy.exp(out[:, 0]), 'r')
plt.plot(ts, numpy.exp(out[:, 1]), 'b')
plt.plot(ts, numpy.exp(out[:, 0] - numpy.sqrt(var[:, 0])), 'r--')
plt.plot(ts, numpy.exp(out[:, 1] - numpy.sqrt(var[:, 1])), 'b--')
plt.plot(ts, numpy.exp(out[:, 0] + numpy.sqrt(var[:, 0])), 'r--')
plt.plot(ts, numpy.exp(out[:, 1] + numpy.sqrt(var[:, 1])), 'b--')
plt.legend(['rabbits', 'wolves'])
plt.gcf().set_size_inches((10, 6))
plt.show()

#%%
out2 = hd.approx(1.5)

plt.plot(ts, out[:, 0], 'g')
plt.plot(ts, out[:, 1], 'y')
plt.plot(ts, out2[:, 0], 'r--')
plt.plot(ts, out2[:, 1], 'b--')
plt.legend(['rabbits', 'wolves'])
plt.show()
#%%
