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

def makeK(x1_, x2_):
    x1, x2 = numpy.meshgrid(x1_, x2_)
    #k = numpy.exp(-0.5 * (x1.T - x2.T)**2 / 4.0)
    k = numpy.min([x1.T, x2.T], axis = 0)
    return k + numpy.eye(len(x1_), len(x2_)) * 1e-10

xs = numpy.linspace(0.0, 1.0, N)

xd = numpy.array([0.0, 1.0])
f = 10.0 * numpy.array([0.0, 0.0])

#xd = numpy.array([1.0, 2.0, 3.5])
#f = 10.0 * numpy.array([0.0, 1.0, 2.0])

Kss = makeK(xs, xs)
Ksd = makeK(xs, xd)
Kds = Ksd.T
Kdd = makeK(xd, xd)

m = Ksd.dot(numpy.linalg.solve(Kdd, f))
Kt = Kss - Ksd.dot(numpy.linalg.solve(Kdd, Kds))
L = numpy.linalg.cholesky(Kt + numpy.eye(Kt.shape[0], Kt.shape[1]) * 1e-10)

for r in range(5):
    plt.plot(xs, m + L.dot(numpy.random.randn(N)), '-')

plt.plot(xd, f, '*')

plt.show()