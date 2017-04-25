#%%
import numpy
import matplotlib.pyplot as plt

xs = numpy.linspace(0.0, 1.0, 5)
ys = numpy.linspace(0.0, 1.0, 10)

xs, ys = numpy.meshgrid(xs, ys)

plt.imshow(xs)
plt.show()
plt.imshow(ys)
plt.show()
#%%

N = 1000

eta2 = 0.001
l2 = 0.05
sigma2 = 0.01

def makeK(x1_, x2_):
    x1, x2 = numpy.meshgrid(x1_, x2_)
    k = eta2 * numpy.exp(-(x1.T - x2.T)**2 / l2)
    #k = numpy.min([x1.T, x2.T], axis = 0)
    return k# + numpy.eye(len(x1_), len(x2_)) * 1e-16#2

xs = numpy.linspace(0.0, 1.0, N)

#xd = numpy.array([0.0, 1.0])
#f = 10.0 * numpy.array([0.0, 0.0])

xd = numpy.linspace(0.0, 1.0, 10)#numpy.array([1.0, 2.0, 3.5])
f = numpy.sin(2 * numpy.pi * xd)#10.0 * numpy.array([0.0, 1.0, 2.0])

Kss = makeK(xs, xs)
Ksd = makeK(xs, xd)
Kds = Ksd.T
Kdd = makeK(xd, xd) + numpy.eye(len(xd)) * sigma2
m = Ksd.dot(numpy.linalg.solve(Kdd, f))
Kt = Kss - Ksd.dot(numpy.linalg.solve(Kdd, Kds))
L = numpy.linalg.cholesky(Kt + numpy.eye(Kt.shape[0], Kt.shape[1]) * 1e-10)

for r in range(100):
    plt.plot(xs, m + L.dot(numpy.random.randn(N)), 'r-', alpha = 0.1)

plt.plot(xs, m, 'b')

plt.plot(xd, f, '*')
plt.title('$\eta^2 = {0}$, $l^2 = {1}$, $\sigma^2 = {2}$'.format(eta2, l2, sigma2), fontsize = 32)
plt.ylim((-2.0, 2.0))

plt.show()
#%%

N = 1000

eta2 = 1.0
l2 = 0.05
sigma2 = 0.00

def makeK(x1_, x2_):
    x1, x2 = numpy.meshgrid(x1_, x2_)
    k = eta2 * numpy.exp(-(x1.T - x2.T)**2 / l2)
    #k = numpy.min([x1.T, x2.T], axis = 0)
    return k# + numpy.eye(len(x1_), len(x2_)) * 1e-16#2

xs = numpy.linspace(0.0, 1.0, N)

#xd = numpy.array([0.0, 1.0])
#f = 10.0 * numpy.array([0.0, 0.0])

xd = numpy.array([0.025, 0.125, 0.3, 0.95])
f = numpy.array([0.0, -1.5, 0.5, 0.0])

Kss = makeK(xs, xs)
Ksd = makeK(xs, xd)
Kds = Ksd.T
Kdd = makeK(xd, xd) + numpy.eye(len(xd)) * sigma2
m = Ksd.dot(numpy.linalg.solve(Kdd, f))
Kt = Kss - Ksd.dot(numpy.linalg.solve(Kdd, Kds))
L = numpy.linalg.cholesky(Kt + numpy.eye(Kt.shape[0], Kt.shape[1]) * 1e-10)

for r in range(3):
    plt.plot(xs, m + L.dot(numpy.random.randn(N)))#, 'r-')

#plt.plot(xs, m, 'b')

plt.plot(xd, f, '*', markersize = 25)
plt.ylim((-2.0, 2.0))

plt.show()