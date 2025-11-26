%% This is a simple example of reconstructing the state variables when the output is already 
% available. Function smd() and obs() are embedded in this file itself. 

clc; 
clear; close all;
%% System properties
m=10; % kg
k=1000; % N/m
c=10; % Ns/m
u=0;
%% System definition
A =[1 0;...
    0 1]; % Fill in 
B = 1; % Fill in
C = [1 0];
p_o=[-2.4, -2.5]; % fill in: Observer pole locations
Ke=; % Ke= transpose of the controller ...
                        % gain matrix of the dual system...
                        % can use place() command 

t_sim=0:.01:10;

%% System calculations
x0=[-2,0]; % released from x=-2 at rest
[t,x]=ode45(@(t,x)smd(x,A,B,u),t_sim,x0);
y=(C*x')'; %([1x2]*[mx2]')' for m time instances
%% Observer calcualtions
x0=[-1,0];
[t_o,x_t] = ode45(@(t_o,x_t)obs(t_o,x_t,A,B,C,u,Ke,y,t),t_sim,x0);


%% Plotting section
% plot the results of x_t and x on the same plots and you can see 
% the performance of the observer. Change the pole locations of the 
% observer and see the difference in behavior of the estimated variables.



%% System function
% spring mass damper state equations
function dx=smd(x,A,B,u)
% state equations
dx = A*x+B*u;
end

%% observer func
function dx_t=obs()
y=interp1(t,Y,t_o); % use this function to find the value of Y as a functoin of t at the instant t_o ...
			% You can also print t_o and see how ODE45 behaves when performing the integration 


dx_t=; % fill in the x_tilda equation here 
end
