# poisson2d 

This is a simple framework for solving PDE on 2D regular grid.

The given example gives the numerical solution to a 2D Laplace equation:

```latex
\Delta f = 0, f: [0, a]\times[0, b]\rightarrow R
```

with the dirichlet boundary conditions:

    f = 0,  if x = 0, y = 0, x = a;
    f = w_0\sin\frac{\pi x}{a}, if y = b

The known analytic solution to this problem is

```latex
f(x,y) = \frac{w_0}{\sinh\frac{\pi b}{a}}\sin\frac{\pi x}{a}\sinh\frac{\pi y}{a}
```
And the result is well fitted to the analytic one.

![image] (comparison.png)

## Heat Equation

Add a new example of solving heat equation on 2D rectangular domain.

```latex
\frac{\partial u}{\partial t}=D\Delta u
```

![image] (diffusion.gif)