function [x,v,fmin,equ]=KKT(Q,c,A,b)
n=size(Q,1);
d=size(A,1);
AA=zeros(n+d,n+d);
bb=zeros(n+d,1);
AA(1:n,1:n)=Q;
AA(n+1:end,1:n)=A;
AA(1:n,n+1:end)=A';
bb(1:n)=-c;
bb(n+1:end)=b;
% disp(AA)
% disp(bb)
res=AA\bb;
disp(res)
x=res(1:n);
v=res(n+1:end);
fmin=1/2*x'*Q*x+c'*x;
equ=A*x-b;