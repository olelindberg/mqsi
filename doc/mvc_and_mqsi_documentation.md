## Mathematical describtion of curves

The curve as a function of arc-length is 
$$
\mathbf{c}(s) = [x(s),y(s)]^T, \quad s \in [a,b]
$$
tangential vector
$$
\mathbf{t}(s) = \mathbf{c}'(s) = [x'(s),y'(s)]^T
$$
arch length parameter
$$
s=\int_a^t \sqrt{x'(t)^2 + y'(t)^2} ~dt

$$

curvature
$$
\kappa = ||\mathbf{t}'(s)|| = ||\mathbf{c}''(s)|| = \sqrt{x''(s)^2 + y''(s)^2}  
$$

## Minimal curvature and minimal curvature variation curves

curvature minimization
$$
\min\int_a^b \kappa^2 ~ds
$$
curvature change minimization
$$
\min\int_a^b \left(\frac{d \kappa}{d s} \right)^2 ~ds
$$

where
$$
\frac{d \kappa }{ds} = \frac{x''(s)\,x'''(s) + y''(s)\,y'''(s)}{\sqrt{{x''(s)}^2 + {y''(s)}^2}}
$$


## Curve discretization

$$
 \{s_i \in [a,b] ~|~ i=0,1,\ldots,N ~ \land ~  s_{i+1}> s_i ~ \}
$$

## Numerical curve representation


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

The curve coordinates are represented by Hermite polynomials
$$
\begin{align}
x(s) &\approx x_{h,i}(s) = H_{x,i}(s), \nonumber  \\
y(s) &\approx y_{h,i}(s) = H_{y,i}(s), \quad s \in [s_i,s_{i+1}], \quad i=1,\ldots,N-1
\end{align}
$$
and the numerical curve representation is
$$
\mathbf{c}(s) \approx \mathbf{c}_{h,i}(s) = [x_{h,i}(s),y_{h,i}(s)]^T, \quad s \in [s_i,s_{i+1}], \quad i=1,\ldots,N -1
$$

## Curve reparametrization

Coordinate transforms
$$
s(t) = \Delta s_i t + s_i 
\quad \Leftrightarrow \quad
t(s) = \frac{s - s_i}{ \Delta s_i} 
$$

Derivatives of the coordinate transforms
$$
\frac{d s}{d t} = \Delta s_i
$$

$$
\frac{d t}{d s} = \frac{1}{\Delta s_i}
$$

$$
\frac{d^2 s}{d t^2} = 0
$$
$$
\frac{d^2 t}{d s^2} = 0
$$

Curve parameter derivatives

$$
\frac{d x}{d s} = \frac{d t}{d s} \frac{d x}{d t} = \frac{1}{\Delta s_i} \frac{d x}{d t}
$$

$$
\frac{d^2 x}{d s^2} = \frac{d}{d s}\left(\frac{d t}{d s} \frac{d x}{d t}  \right) = \frac{d^2 t}{d s^2} \frac{d x}{d t} + 
\left(\frac{d s}{d t} \right)^2 \frac{d^2 x}{d s^2}
=
\left(\frac{1}{\Delta s_i} \right)^2 \frac{d^2 x}{d s^2}
$$


## Minimization solver 
###  Gradient decent
$$
\mathbf{q}^{n+1} = \mathbf{q}^n - \frac{d f}{d \mathbf{a}}
$$

In the curvature variation integral the integrad is
$$
f(s,a) = 
\left( 
\frac{x''(s,a)\, x'''(s,a) + y''(s,a)\, y'''(s,a)}{\sqrt{x''(s,a)^2 + y''(s,a)^2}} 
\right) \ 2
$$

T \he partial derivative of the integrad with respect to the parameter $ a \in \{f_i,f'_i,f''_i,f_{i+1},f'_{i+1},f''_{i+1} \}$

  $$
\frac{\partial f}{\partial a} 
=   
\frac{2u}{v^3} \left( v \ u_a - u \ v_a \right)
$$

where 

$$
\begin{align}
u   &= x'' x''' + y'' y''' \\
u_a &= x''_a x''' + x'' x'''_a + y''_a y''' + y'' y'''_a \\
v   &= \sqrt{x''^2 + y''^2} \\
v_a &= \frac{x'' x''_a + y'' y''_a}{\sqrt{x''^2 + y''^2}}
\end{align}
$$
## Curve parametrization
$$
s(\xi) = \frac{\xi + 1}{2}(s_{i+1} - s_i)
$$

$$
\frac{d s}{d \xi} = \frac{1}{2}(s_{i+1} - s_i)
$$

$$
\frac{d x}{d \xi} = \frac{d s}{d \xi} \frac{d x}{d s}
$$

$$
\frac{d^2 x}{d \xi^2} = \frac{d}{d \xi}\left(\frac{d s}{d \xi} \frac{d x}{d s}  \right)
$$

$$
\frac{d^2 x}{d \xi^2} = 
\frac{d^2  s}{d \xi^2} \frac{d x}{d s}
+
\left(\frac{d s}{d \xi} \right)^2 \frac{d^2 x}{d s^2}
$$

## Coordinate transforms

$$
s = s(t), \quad t \in [a,b]
$$

$$
\begin{align}
\frac{d x}{d s}     &= \frac{d t}{d s} \frac{d x}{d t} \\
\frac{d^2 x}{d s^2} &= \frac{d^2 t}{d s^2} \frac{d x}{d t} + \left(\frac{d t}{d s} \right)^2 \frac{d^2 x}{d t^2} \\
\frac{d^3 x}{d s^3} &= \frac{d^3  t}{d s^3} \frac{d x}{d t} + 3\frac{d^2  t}{d s^2} \frac{d  t}{d s} \frac{d^2 x}{d t^2} + \left(\frac{d t}{d s} \right)^3 \frac{d^3 x}{d t^3} \\
\end{align}


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
