function f=Rosenbrock(x)
x_odd=x(1:2:end);
x_even=x(2:2:end);
f=sum(100*(x_odd.^2-x_even).^2+(x_odd-1).^2);