# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 10:50:22 2018

@author: Masha
"""

#однопараметрический анализ по параметру k1
import matplotlib.pyplot as plot
import math


_k1 = 0.12
_k1m = 0.005
_k1m1 = [0.001, 0.005, 0.01, 0.015, 0.02]
_k3 = 0.0032
_k2 = 1.05
_k3m = 0.002
_k3m1 =  [0.0005, 0.001, 0.002, 0.003, 0.004]
    
alkp = math.sqrt(_k3) / (math.sqrt(_k3m) + math.sqrt(_k3))
x = [i / 1000.0 for i in range(987)]
n2 = len(x)
y = [(1 - i) * alkp for i in x]
z = [1 - xx - yy for xx, yy in zip(x, y)]

for k1m in _k1m1:
    k1 = []
    sp = []
    det = []
    DI = []

    xh = []
    yh = []
    k1h = []

    xsn = []
    ysn = []
    k1sn = []

    xdi = []
    ydi = []
    k1di = []


    k1.append((k1m*x[0] + _k2*x[0]*z[0]**2) / z[0])
    a11 = -k1[0] - k1m - _k2*z[0]**2 + 2*_k2*x[0]*z[0]
    a12 = -k1[0] + 2*_k2*x[0]*z[0]
    a21 = -2* _k3 * z[0]
    a22 = -2*_k3 * z[0] - 2*_k3m*y[0]

    sp.append(a11 + a22)
    det.append(a11*a22 - a12*a21)
    DI.append(sp[0]**2 - 4*det[0])

    for i in range(1, n2):
        k1.append((k1m*x[i] + _k2*x[i]*z[i]**2)/z[i])
        a11 = -k1[i] - k1m - _k2*z[i]**2 + 2*_k2*x[i]*z[i]
        a12 = -k1[i] + 2*_k2*x[i]*z[i]
        a21 = -2* _k3 * z[i]
        a22 = -2*_k3 * z[i] - 2*_k3m*y[i]
        sp.append(a11 + a22)
        det.append(a11*a22 - a12*a21)
        DI.append(sp[i]**2 - 4*det[i])
    
        if sp[i]*sp[i-1] <= 0:
            yh.append(y[i])
            xh.append(x[i])
            k1h.append(k1[i])

        if det[i]*det[i-1] <= 0:
            ysn.append(y[i])
            xsn.append(x[i])
            k1sn.append(k1[i])

        if DI[i]*DI[i-1] <= 0:
            ydi.append(y[i])
            xdi.append(x[i])
            k1di.append(k1[i])


    plot.figure(2)
    plot.subplot(111)
    plot.title('Parametric Portrait')
    plot.xlabel('k1')
    plot.ylabel('y1, y2')
    plot.plot(k1, x, 'g-', label = 'y1')
    plot.plot(k1, y, 'r-', label = 'y2')
    plot.plot(k1h, xh, 'g*', label = 'Hopf')
    plot.plot(k1h, yh, 'r*', label = 'Hopf')
    plot.plot(k1sn, xsn, 'g+', label = 'Saddle-node')
    plot.plot(k1sn, ysn, 'r+', label = 'Saddle-node')
    plot.legend()
    plot.text(0.09, 1.01, 'k values')
    plot.grid()
    plot.title('Dependence on Parameter k₋₁ = '+ str(k1m))
    plot.axis([0.05, 0.3, 0.0, 1.10])
    plot.savefig('k1_' + str(k1m) + '.png')
    plot.show()
    
for k3m in _k3m1:
    alkp = math.sqrt(_k3) / (math.sqrt(k3m) + math.sqrt(_k3))
    y = [(1 - i) * alkp for i in x]
    z = [1 - xx - yy for xx, yy in zip(x, y)]
    k1 = []
    sp = []
    det = []
    DI = []

    xh = []
    yh = []
    k1h = []

    xsn = []
    ysn = []
    k1sn = []

    xdi = []
    ydi = []
    k1di = []


    k1.append((_k1m*x[0] + _k2*x[0]*z[0]**2) / z[0])
    a11 = -k1[0] - _k1m - _k2*z[0]**2 + 2*_k2*x[0]*z[0]
    a12 = -k1[0] + 2*_k2*x[0]*z[0]
    a21 = -2* _k3 * z[0]
    a22 = -2*_k3 * z[0] - 2*k3m*y[0]

    sp.append(a11 + a22)
    det.append(a11*a22 - a12*a21)
    DI.append(sp[0]**2 - 4*det[0])

    for i in range(1, n2):
        k1.append((_k1m*x[i] + _k2*x[i]*z[i]**2)/z[i])
        a11 = -k1[i] - _k1m - _k2*z[i]**2 + 2*_k2*x[i]*z[i]
        a12 = -k1[i] + 2*_k2*x[i]*z[i]
        a21 = -2* _k3 * z[i]
        a22 = -2*_k3 * z[i] - 2*k3m*y[i]
        sp.append(a11 + a22)
        det.append(a11*a22 - a12*a21)
        DI.append(sp[i]**2 - 4*det[i])
    
        if sp[i]*sp[i-1] <= 0:
            yh.append(y[i])
            xh.append(x[i])
            k1h.append(k1[i])

        if det[i]*det[i-1] <= 0:
            ysn.append(y[i])
            xsn.append(x[i])
            k1sn.append(k1[i])

        if DI[i]*DI[i-1] <= 0:
            ydi.append(y[i])
            xdi.append(x[i])
            k1di.append(k1[i])


    plot.figure(2)
    plot.subplot(111)
    plot.title('Parametric Portrait')
    plot.xlabel('k1')
    plot.ylabel('y1, y2')
    plot.plot(k1, x, 'g-', label = 'y1')
    plot.plot(k1, y, 'r-', label = 'y2')
    plot.plot(k1h, xh, 'g*', label = 'Hopf')
    plot.plot(k1h, yh, 'r*', label = 'Hopf')
    plot.plot(k1sn, xsn, 'g+', label = 'Saddle-node')
    plot.plot(k1sn, ysn, 'r+', label = 'Saddle-node')
    plot.legend()
    plot.text(0.09, 1.01, 'k values')
    plot.grid()
    plot.title('Dependence on Parameter k₋3 = '+ str(k3m))
    plot.axis([0.05, 0.2, 0.0, 1.07])
    plot.savefig('k3_' + str(k3m) + '.png')
    plot.show()