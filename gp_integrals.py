#%%

import numpy
from numpy import zeros
from numpy import eye as identity
from math import erf, exp, pi, sqrt
import matplotlib.pyplot as plt

l = 0.5
a = 0.5
sig = 0.0

def R(tj, tk):
    return a * exp(-(tj - tk)**2 / (2 * l**2))

def QR(tj, tk):
    return a * pi * l**2 * (erf((tj - tk) / (2 * l)) + erf(tk / (2 * l)))

def RQ(tj, tk):
    return a * QR(tk, tj)

def QQ(tj, tk):
    t1 = exp(-tj**2 / (4 * l**2)) - exp(-(tj - tk)**2 / (4 * l**2)) + exp(-tk**2 / (4 * l**2))
    t2 = tj * erf(tj / (2 * l))
    t3 = (tk - tj) * erf((tj - tk) / (2 * l))
    t4 = tk * erf(tk / (2 * l))
    return a * pi * l**2 * (2 * l * t1 / sqrt(pi) + t2 + t3 + t4)

N = 10
Nt = 100

ts = [0] * N

I = 20

xs = numpy.linspace(0, N, N)
y = numpy.sin(5 * pi * xs / N)#numpy.random.randn(N)
ys = numpy.linspace(0, N, Nt)

def cov(lx, ly, xs, ys):
    K = zeros((len(xs), len(ys)))

    for i in range(len(xs)):
        for j in range(len(ys)):
            if lx[i] == 0 and lx[j] == 0:
                K[i, j] = R(xs[i], ys[j])
            elif lx[i] == 1 and ly[j] == 0:
                K[i, j] = QR(xs[i], ys[j])
            elif lx[i] == 0 and ly[j] == 1:
                K[i, j] = RQ(xs[i], ys[j])
            else:
                K[i, j] = QQ(xs[i], ys[j])

    return K

xs = numpy.concatenate(([0], xs))
y = numpy.concatenate(([0], y))
lx = [1] + [0] * N
lyd = [0] * Nt
lyi = [1] * Nt

K = cov(lx, lx, xs, xs)
KsK = cov(lyd, lx, ys, xs)
KKs = KsK.T
KsKs = cov(lyd, lyd, ys, ys)

KsKi = cov(lyi, lx, ys, xs)
KKsi = KsKi.T
KsKsi = cov(lyi, lyi, ys, ys)

def mu(K, KsK, y):
    Kt = K.copy()

    for i in range(N):
        Kt[i + 1, i + 1] += sig**2;

    return KsK.dot(numpy.linalg.solve(Kt, y))

def cov(K, KsK, KsKs):
    KKs = KsK.T

    Kt = K.copy()

    for i in range(N):
        Kt[i + 1, i + 1] += sig**2;

    return KsKs - KsK.dot(numpy.linalg.solve(Kt, KKs))

yd = numpy.random.multivariate_normal(mu(K, KsK, y), cov(K, KsK, KsKs), I)
yi = numpy.random.multivariate_normal(mu(K, KsKi, y), cov(K, KsKi, KsKsi), I)

plt.imshow(K)
plt.colorbar()
plt.show()

for i in range(I):
    plt.plot(ys, yd[i], '-*', alpha = 0.25)
plt.plot(ys, mu(K, KsK, y), 'r')
plt.plot(xs[1:], y[1:], 'b*')
plt.show()

for i in range(I):
    plt.plot(ys, yi[i], '*', alpha = 0.25)
plt.plot(ys, mu(K, KsKi, y), 'g')