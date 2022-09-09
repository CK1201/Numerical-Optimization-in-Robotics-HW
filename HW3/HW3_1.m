clc,clear
close all
%%
Q=[1,2;3,4];
c=[5;6];
A=[7,8];
b=[9];
[x,v,fmin,equ]=KKT(Q,c,A,b);