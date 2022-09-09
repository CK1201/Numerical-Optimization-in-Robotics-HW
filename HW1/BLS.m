% Backtracking Line Search
function [xk,f_best,iter_times,x_list,f_list]=BLS(fun,grad,x0,c,tau0)
xk=x0;
iter_times=0;
x_list=xk;
f_list=[];
while(norm(grad(xk))>10e-5)
    iter_times=iter_times+1;
    d=-grad(xk);
    tau=tau0;
    fk=fun(xk);
    f_list=[f_list,fk];
    while(fun(xk+tau*d)>fk+c*tau*d'*grad(xk))
        tau=tau/2;
    end
    xk=xk+tau*d;
    x_list=[x_list;xk];
end
f_best=fun(xk);
f_list=[f_list,f_best];