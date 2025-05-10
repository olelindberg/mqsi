from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt

from hermite_quintic        import hermite_quintic
from curvature              import curvature


curve = "wicket3"


t = np.arange(-1, 1.05, 0.05)
H0  = hermite_quintic(0,t)
Ht  = hermite_quintic(1,t)
Htt = hermite_quintic(2,t)


if curve== "wicket2":

    r     = 1.0
    angle = np.pi/2
    k     = 1/r
    cx    = [r,  0, -r*angle**2, -r,   0, r*angle**2]
    cy    = [0,r*angle, 0, 0,-r*angle, 0]

elif curve == "wicket3":

    r     = 1.0
    angle = np.pi/4
    k     = 1/r
    cx    = [r,       0, -r*angle**2, 0, -r*angle,           0]
    cy    = [0, r*angle,           0, 1,        0, -r*angle**2]


#    denom = (x_t**2 + y_t**2)**1.5
#    return (x_t * y_tt - y_t * x_tt)/denom


xx = r*np.cos(angle*(t+1))
yy = r*np.sin(angle*(t+1))

xx_t = -r*angle*np.sin(angle*(t+1))
yy_t = r*angle*np.cos(angle*(t+1))

xx_tt = -r*angle**2*np.cos(angle*(t+1))
yy_tt = -r*angle**2*np.sin(angle*(t+1))


x  = H0@cx
y  = H0@cy

x_t  = Ht@cx
y_t  = Ht@cy

x_tt  = Htt@cx
y_tt  = Htt@cy

k  = curvature(x_t, y_t, x_tt,y_tt)
kk = curvature(xx_t, yy_t, xx_tt,yy_tt)

print("radius = ",r)

fig, axs = plt.subplots(2,3, figsize=(10, 8))
axs[0,0].plot(t, x)
axs[0,0].plot(t, xx)
axs[0,0].set_title('x')
axs[0,0].grid(True)

axs[1,0].plot(t, y)
axs[1,0].plot(t, yy)
axs[1,0].set_title('y')
axs[1,0].grid(True)

axs[0,1].plot(t, x_t)
axs[0,1].plot(t, xx_t)
axs[0,1].set_title('x_t')
axs[0,1].grid(True)

axs[1,1].plot(t, y_t)
axs[1,1].plot(t, yy_t)
axs[1,1].set_title('y_t')
axs[1,1].grid(True)

axs[0,2].plot(t, x_tt)
axs[0,2].plot(t, xx_tt)
axs[0,2].set_title('x_tt')
axs[0,2].grid(True)

axs[1,2].plot(t, y_tt)
axs[1,2].plot(t, yy_tt)
axs[1,2].set_title('y_tt')
axs[1,2].grid(True)

plt.figure()
plt.plot(t, k)
plt.plot(t, kk)
plt.title('curvature')
plt.grid(True)

plt.figure()

plt.plot(xx, yy,'k.-')
plt.plot(x, y,'r.-')
plt.title('x-y')
plt.grid(True)
plt.axis('equal')

plt.show()

