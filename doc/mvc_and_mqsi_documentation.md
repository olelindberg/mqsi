## Mathematical describtion of curves

The curve as a function of arc-length is 
$$
\mathbf{c}(s) = [x(s),y(s)]^T
$$
tangential vector
$$
\mathbf{t}(s) = \mathbf{c}'(s) = [x'(s),y'(s)]^T
$$
arch length 
$$
l=\int_a^b \sqrt{x'(s)^2 + y'(s)^2} ~ds

$$

curvature
$$
\kappa = ||\mathbf{t}'(s)|| = ||\mathbf{c}''(s)|| = \sqrt{x''(s)^2 + y''(s)^2}  
$$
curvature minimization
$$
\min\int_a^b \kappa^2 ~ds
$$
curvature change minimization
$$
\min\int_a^b \left(\frac{d \kappa}{d s} \right)^2 ~ds
$$

## Numerical approximation

Curve discretization

s 


## Minimization solver 
###  Gradient decent
$$
\mathbf{q}^{n+1} = \mathbf{q}^n - \frac{d f}{d \mathbf{a}}
$$


## Curve parametrization
$$
s(\xi) = \frac{\xi + 1}{2}(s_{i+1} - s_i)
$$

$$
\frac{\partial{s}}{\partial \xi} = \frac{1}{2}(s_{i+1} - s_i)
$$


$$
\frac{\partial{x}}{\partial \xi} = \frac{\partial{s}}{\partial \xi} \frac{\partial{x}}{\partial s}
$$

$$
\frac{\partial^2 x}{\partial \xi^2} = \frac{\partial}{\partial \xi}\left(\frac{\partial{s}}{\partial \xi} \frac{\partial{x}}{\partial s}  \right)
$$

$$
\frac{\partial^2 x}{\partial \xi^2} = 
\frac{\partial^2 {s}}{\partial \xi^2} \frac{\partial x}{\partial s}
+
\left(\frac{\partial{s}}{\partial \xi} \right)^2 \frac{\partial^2 x}{\partial s^2}
$$



## Monotone quintic spline interpolation

The theory for monotone quintic spline interpolation is described in 

Walter Hess, Jochen W. Schmidt*, "Positive quartic, monotone quintic C2-spline interpolation in one and two dimensions", Journal of Computational and Applied Mathematics 55 (1994) 51-67.

For s to be monotone on I, are equivalent to
$$
\begin{align}
p_i                                          & \geq 0, \quad i = 0,\ldots,n \\ 
4 p_{i-1} + h_i P_{i-1}                      & \geq 0, \quad i = 1,\ldots,n\\ 
4 p_i - h_i P_i                                   & \geq 0, \quad i = 1,\ldots,n\\
60 (z_i- z_{i-1}) - 24 h_i(p_i +p_{i-1}) + 3h_i^2(P_i - P_{i-1}) & \geq 0, \quad i = 1,\ldots,n
\end{align}
$$

# Applications


## 2D Wicket
The two-dimensional wicket is a half circle
$$
\begin{align}
x &= r \cos(\theta), \\ 
y &= r \sin(\theta),\quad \theta \in [0,\pi].
\end{align}
$$
Since the MVC equations are formulated in terms of arc length $s$ we are using the relation between angle and arc length of a circle 
$$
\theta = \frac{s}{r}
$$
to express the wicket in terms of arc length:
$$
\begin{align}
x &= r \cos \left( \frac{s}{r} \right), \\ 
y &= r \sin \left( \frac{s}{r} \right),\quad s \in [0,\pi \, r].
\end{align}
$$
For the inital condition we need the derivatives with respect to arc length $x'$, $y'$, $x''$ and $y''$. The first derivatives are
$$
\begin{align}
x' &= \cos \left( \frac{s}{r} \right), \\ 
y' &= \sin \left( \frac{s}{r} \right),\quad s \in [0,\pi \, r].
\end{align}
$$
and the second derivatives are
$$
\begin{align}
x'' &= \cos \left( \frac{s}{r} \right), \\ 
y'' &= \sin \left( \frac{s}{r} \right),\quad s \in [0,\pi \, r].
\end{align}
$$

