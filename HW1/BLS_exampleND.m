clc,clear
close all
%% call
N=10;
% x0=rand(1,N)*2*2-2;
x0=ones(1,N);
x0(1:2:end)=-1;
x0=x0*2;
tao_list=0.01:0.1:1;
c_list=[1e-5,1e-4,1e-3,1e-2,1e-1];
exp_size=length(c_list);
% exp_size=length(tao_list);
csm_time=zeros(1,exp_size);
iter=zeros(1,exp_size);
for i=1:exp_size
%     tic
    [x_best,f_best,iter(i)]=BLS(@Rosenbrock,@Rosenbrock_grad,x0,c_list(i),1);
%     csm_time(i)=toc;
end

% plot(tao_list,iter,'.-')
semilogx(c_list,iter,'.-')
grid on
ylabel('迭代次数')
xlabel('常数c')