m=1;
M=5;
L=2;
g=-9.81;
up=1;
d=15;
%ref=[1;0;pi;0];
t_sim=0:0.001:30;x0=[0;0;pi+0.1;0];
[t,x]=ode45(@(t,x)CtP(x,M,m,L,g,0,d),t_sim,x0);
plot(t,rad2deg(x(:,3)))
% CtP_render