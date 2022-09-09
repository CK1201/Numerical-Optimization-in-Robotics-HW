# Backtracking Line Search

An unconstrained optimizer.

## How to use

It can be used to optimize any smooth objective function as long as the closed form of the function and the gradient of the function is given.

input: objective function, gradient of the objective function, initial guess, constant c, initial tau

output: best x, best fun value, iteration times, history x, history fun value

```matlab
[xk,f_best,iter_times,x_list,f_list]=BLS(fun,grad,x0,c,tau0)
```

See "BLS_example2D.m" and “BLS_exampleND.m” for the calling procedure. 
