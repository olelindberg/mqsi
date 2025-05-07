## Mathematical describtion of curves

The curve as a function of arc-length is 
$$
\mathbf{c}(s) = [x(s),y(s)]^T, \quad s \in [a,b]
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

## Curve discretization

$$
 \{s_i \in [a,b] ~|~ i=0,1,\ldots,N ~ \land ~  s_{i+1}> s_i ~ \}
$$

## Solution approximation


Hermite quintic basis functions with $t \in [0,1]$
$$
\begin{align}
h_0​(t)&=1−10t^3+15t^4−6t^5     \nonumber \\
h_1(t)&=t−6t^3+8t^4−3t^5       \nonumber \\
h_2(t)&=\tfrac{1}{2}t^2−\tfrac{3}{2}t^3+\tfrac{3}{2}t^4−\tfrac{1}{2}t^5 \nonumber \\
h_3(t)&=10t^3−15t^4+6t^5       \nonumber \\
h_4(t)&=−4t^3+7t^4−3t^5        \nonumber \\
h_5(t)&=\tfrac{1}{2}t^3−t^4+\tfrac{1}{2}t^5
\end{align} 
$$

Parameter
$$
t = \frac{s-s_i}{\Delta s_i}, \quad s \in [s_i,s_{i+1}]
$$

Hermite quintic polynomial
$$
\begin{align}
H_i(s)=&f_ih_0​(t)+\Delta s_i f_i'h_1​(t)+\Delta s_i ^2f_i''h_2​(t) \nonumber \\
+&f_{i+1}h_3​(t)+\Delta s_i f_{i+1}'h_4​(t)+\Delta s_i ^2f_{i+1}''h_5​(t), \quad s \in [s_i,s_{i+1}]
\end{align}
$$


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