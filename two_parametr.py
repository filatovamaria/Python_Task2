import numpy
import matplotlib.pyplot as plot
import math
from sympy import lambdify, Symbol
from scipy.integrate import odeint

# K1K2
def K1_sp(x, y, z, k2, k3, k3m):    
    return 2 * (k2 * x**2 * z - k3 * z - k3m * y) / (x + z)

def K1_det(x, y, z, k2, k3, k3m):
    return 2 * k2 * x * z * (k3*z*(x - 1) + k3m * x * y) / (k3*z*(x + z - 1) + k3m * y *(x + z))

def k1m_(x, z, k1, k2):
    return z*(k1 - k2 * x * z) / x


_k1 = 0.12
_k1m= 0.005
_k3 = 0.0032
_k2 = 1.05
_k3m= 0.002


alkp = math.sqrt(_k3) / (math.sqrt(_k3m) + math.sqrt(_k3))

x = [i / 1000.0 for i in range(1, 999)]
y = [(1 - _x) * alkp for _x in x]
z = [1 - _y - _x for _y, _x in zip(y, x)]

# расчет линии кратности
K1 = [K1_sp(_x, _y, _z, _k2, _k3, _k3m) for _x, _y, _z in zip(x, y, z)]
K1m = [k1m_(_x, _z, _k1, _k2) for _x, _z, _k1 in zip(x, z, K1)]

# расчет линии нейтральности
k1 = [K1_det(_x, _y, _z, _k2, _k3, _k3m) for _x, _y, _z in zip(x, y, z)]
k1m = [k1m_(_x, _z, _k1, _k2) for _x, _z, _k1 in zip(x, z, k1)]


# построение параметрического портрета
plot.figure(1)
plot.subplot(111)
plot.title('Parametric Portrait')
plot.xlabel('k1m')
plot.ylabel('k1')
plot.plot(k1m, k1, 'g.', label = 'Saddle-node')
plot.plot(K1m, K1, 'r-', label = 'Hopf')
plot.legend()
plot.grid()
plot.title('Parametric Portrait')
plot.axis([-0.01, 0.022, 0, 0.22])

plot.annotate (u'I', xy=(-0.002, 0.142), fontsize=18)
plot.annotate (u'I', xy=(0.012, 0.013), fontsize=18)
plot.annotate(u'II', xy=(-0.009, 0.0755), fontsize=14)
plot.annotate(u'II+', xy=(-0.003, 0.06), fontsize=16)
plot.annotate (u'III - oscillation', xy=(0.011, 0.14), fontsize=16)
plot.annotate(u'sn2',xy=(0.0015, 0.055))
plot.annotate(u'sn1',xy=(-0.008, 0.06))
plot.annotate(u'C', xy=(0.0038, 0.097), fontsize=14)
plot.annotate(u'TB', xy=(0.0018, 0.11), fontsize=12)
plot.savefig('param.png')
plot.show()



# построение фазового портрета

k1 = Symbol('k1', Positive=True)
k1m = Symbol('km1', Positive=True)
k2 = Symbol('k2', Positive=True)
k3 = Symbol('k3_0', Positive=True)
k3m = Symbol('k3_0', Positive=True)
x = Symbol("x", Positive=True)
y = Symbol("y", Positive=True)

dxdt = k1*(1 - x - y) - k1m*x - k2*x*(1 - x - y)**2
dydt = k3*(1 - x - y)**2 - k3m* y**2

def dzdt(z, t, k1, k1m, k2, k3, k3m):
    x = z[0]
    y = z[1]
    dxdt = k1*(1 - x - y) - k1m*x - k2*x*(1 - x - y)**2
    dydt = k3*(1 - x - y)**2 - k3m* y**2
    dzdt = [dxdt, dydt]
    return dzdt


t = numpy.linspace(0,50000,100000) 

z1 = [0, 0.00]
res1 = odeint(dzdt, z1, t, args = (_k1,_k1m,_k2,_k3,_k3m))
res1_x = [r[0] for r in res1]
res1_y = [r[1] for r in res1]

z2 = [0.1, 0.3]
res2 = odeint(dzdt, z2, t, args = (_k1,_k1m,_k2,_k3,_k3m))
res2_x = [r[0] for r in res2]
res2_y = [r[1] for r in res2]

z3 = [3, 0.00]
res3 = odeint(dzdt, z3, t, args = (_k1,_k1m,_k2,_k3,_k3m))
res3_x = [r[0] for r in res3]
res3_y = [r[1] for r in res3]

z4 = [0, 3.00]
res4 = odeint(dzdt, z4, t, args = (_k1,_k1m,_k2,_k3,_k3m))
res4_x = [r[0] for r in res4]
res4_y = [r[1] for r in res4]

z5 = [0.015, 0.15]
res5 = odeint(dzdt, z5, t, args = (_k1,_k1m,_k2,_k3,_k3m))
res5_x = [r[0] for r in res5]
res5_y = [r[1] for r in res5]

plot.figure(1)
plot.xlabel('x')
plot.ylabel('y')

plot.plot(res1_x, res1_y, 'r-', label = u'y(x)')
plot.plot(res2_x, res2_y, 'r-')
plot.plot(res3_x, res3_y, 'r-')
plot.plot(res4_x, res4_y, 'r-')
plot.plot(res5_x, res5_y, 'r-')

f1 = lambdify((x, y, k1, k1m, k2), dxdt)
f2 = lambdify((x, y, k3, k3m), dydt)
Y, X = numpy.mgrid[0:.5:100j, 0:1:200j]
U = f1(X, Y, _k1, _k1m, _k2)
V = f2(X, Y, _k3, _k3m)
speed = numpy.sqrt(X**2 + Y**2)
plot.streamplot(X, Y, U, V, density = [2, 0.5], color = 'y')

plot.axis([0, 1, 0, 0.5])
plot.legend()
plot.grid(True)
plot.title('Phase Portrait')
plot.savefig('phase.png')
plot.show()

plot.figure(1)
plot.xlabel('t')
plot.ylabel('x, y')
plot.plot(t, res_x, 'r-', label = u'x(t)')
plot.plot(t, res_y, 'g-', label = u'y(t)')
plot.axis([0, 50000, 0, 1.2])
plot.legend()
plot.grid(True)
plot.title('Dependence on t')
plot.savefig('t_dependence.png')
plot.show()


