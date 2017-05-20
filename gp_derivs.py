#%%

import numpy
from numpy import zeros, linspace, eye
from math import erf, exp, pi, sqrt, sin
from scipy.integrate import odeint
import matplotlib.pyplot as plt

l = 1.0
a = 1.0
s = 1e-1

gl = 9.8

def QQ(tj, tk):
    return exp(-((tj - tk)**2/(2 * l**2))) * a**2

def QR(tj, tk):
    return (a**2 * exp(-((tj - tk)**2/(2 * l**2))) * (tj - tk))/l**2

def RQ(tj, tk):
    return QR(tk, tj)

def RR(tj, tk):
    return (a**2 * exp(-((tj - tk)**2/(2 * l**2))))/l**2 - (a**2 * exp(-((tj - tk)**2/(2 * l**2))) * (tj - tk)**2)/l**4

def QT(tj, tk):
    return -((a**2 * exp(-((tj - tk)**2/(2 * l**2))))/l**2) + (a**2 * exp(-((tj - tk)**2/(2 * l**2))) * (tj - tk)**2)/l**4

def TQ(tj, tk):
    return QT(tk, tj)

def RT(tj, tk):
    return (3 * a**2 * exp(-((tj - tk)**2/(2 * l**2))) * (tj - tk))/l**4 - (a**2 * exp(-((tj - tk)**2/(2 * l**2))) * (tj - tk)**3)/l**6

def TR(tj, tk):
    return RT(tk, tj)

def TT(tj, tk):
    return (3 * a**2 * exp(-((tj - tk)**2/(2 * l**2))))/l**4 - (6 * a**2 * exp(-((tj - tk)**2/(2 * l**2))) * (tj - tk)**2)/l**6 + (a**2 * exp(-((tj - tk)**2/(2 * l**2))) * (tj - tk)**4)/l**8

def full(y, t0):
    f = numpy.zeros(2)

    f[0] = y[1]
    f[1] = -gl * sin(y[0])

    return f

def linearized(y, t0):
    f = numpy.zeros(2)

    f[0] = y[1]
    f[1] = -gl * y[0]

    return f

T = 5.0
N = 25
Nt = 25
y0 = numpy.zeros(2)
y0[0] = numpy.pi / 4.0
ts = linspace(0, T, N)
tts = linspace(0, T, Nt)
y = odeint(linearized, y0, ts)[:, 0]
plt.plot(ts, y)
y = odeint(full, y0, ts)[:, 0]
plt.plot(ts, y)
plt.legend(['linearized', 'full'])
plt.show()

#%%

I = 15

def cov(f, xs, ys):
    K = zeros((len(xs), len(ys)))

    for i in range(len(xs)):
        for j in range(len(ys)):
            K[i, j] = f(xs[i], ys[j])

    return K

lyi = [0] * N
lydd = [1] * N

K = cov(QQ, ts, ts)
KsK = cov(QQ, tts, ts)
KKs = KsK.T
KsKs = cov(QQ, tts, tts)

KsKi = cov(TQ, tts, ts)
KKsi = KsKi.T
KsKsi = cov(TT, tts, tts)

def mu(K, KsK, y):
    Kt = K.copy()

    for i in range(K.shape[0]):
        Kt[i, i] += s**2;

    return KsK.dot(numpy.linalg.solve(Kt, y))

def cov(K, KsK, KsKs):
    KKs = KsK.T

    Kt = K.copy()

    for i in range(K.shape[0]):
        Kt[i, i] += s**2;

    return KsKs - KsK.dot(numpy.linalg.solve(Kt, KKs))

yd = numpy.random.multivariate_normal(mu(K, KsK, y), cov(K, KsK, KsKs), I)
yi = numpy.random.multivariate_normal(mu(K, KsKi, y), cov(K, KsKi, KsKsi), I)

for i in range(I):
    plt.plot(tts, yd[i], '-*', alpha = 0.25)
plt.plot(tts, mu(K, KsK, y), 'r')
plt.plot(ts, y, 'b*')
plt.show()

for i in range(I):
    plt.plot(tts, yi[i], '*', alpha = 0.25)
plt.plot(tts, mu(K, KsKi, y), 'g')