from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt

from mvc_integrand_jacobian import mvc_integrand_jacobian
from mvc_objective_function import mvc_objective_function
from hermite_quintic        import hermite_quintic
from curvature              import curvature

t = np.arange(-1, 1.05, 0.05)
H0  = hermite_quintic(0,t)
Ht  = hermite_quintic(1,t)
Htt = hermite_quintic(2,t)

r = 4.0
l = np.pi*r
k = 1/r

#    denom = (x_t**2 + y_t**2)**1.5
#    return (x_t * y_tt - y_t * x_tt)/denom


xx = r*np.cos(np.pi/2*(t+1))
yy = r*np.sin(np.pi/2*(t+1))

xx_t = -r*np.pi/2*np.sin(np.pi/2*(t+1))
yy_t = r*np.pi/2*np.cos(np.pi/2*(t+1))

xx_tt = -r*np.pi/2*np.pi/2*np.cos(np.pi/2*(t+1))
yy_tt = -r*np.pi/2*np.pi/2*np.sin(np.pi/2*(t+1))

cx = [r,  0, -r/4*np.pi**2, -r,   0, r/4*np.pi**2]
cy = [0,r*np.pi/2, 0, 0,-r*np.pi/2, 0]

x  = H0@cx
y  = H0@cy

x_t  = Ht@cx
y_t  = Ht@cy

x_tt  = Htt@cx
y_tt  = Htt@cy

k  = curvature(x_t, y_t, x_tt,y_tt)
kk = curvature(xx_t, yy_t, xx_tt,yy_tt)

print("radius = ",r)
print("length = ",l)

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

