function grad=Rosenbrock_grad(x)
x_odd=x(1:2:end);
x_even=x(2:2:end);
grad=x;
grad(1:2:end)=200*(x_odd.^2-x_even)*2.*x_odd+2*(x_odd-1);
grad(2:2:end)=-200*(x_odd.^2-x_even);