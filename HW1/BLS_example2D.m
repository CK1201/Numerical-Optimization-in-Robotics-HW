clc,clear
close all
%% call
N=2;
x0_range=2;
x0=rand(1,N)*2*x0_range-x0_range;
x0=[-1,1]*x0_range;
tic
[x_best,f_best,iter_times,x_list,f_list]=BLS(@Rosenbrock,@Rosenbrock_grad,x0,1e-3,1);
toc
[x_best1,f_best1,iter_times1,x_list1,f_list1]=BLS(@Rosenbrock,@Rosenbrock_grad,x0,1e-3,0.1);

%% draw
size_x=2.1;
x1=-size_x:0.01:size_x-0.001;
x2=-size_x+2:0.01:size_x-0.001+0.5;
[x1,x2] = meshgrid(x1,x2);
f=zeros(size(x1));
for i=1:size(x1,1)      
    for j=1:size(x1,2)
        f(i,j)=Rosenbrock([x1(i,j),x2(i,j)]);
    end
end
% contour(x1,x2,f)
% mesh(x1,x2,f)
% hold on
% plot(x_list(:,1),x_list(:,2),'.-','MarkerSize',12)
% plot3(1,1,0,'g.','MarkerSize',20)
% for i=1:size(x_list,1)
%     plot3(x_list(i,1),x_list(i,2),f_list(i),'r.','MarkerSize',20)
% end
% plot(x_best(1),x_best(2),'r.','MarkerSize',16)

semilogy(1:iter_times+1,f_list)
grid on
hold on
semilogy(1:iter_times1+1,f_list1)
legend('tau0=1','tau0=0.1')
xlabel('迭代次数')
ylabel('f-f*')