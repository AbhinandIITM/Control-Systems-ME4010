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
A =[0 1;...
    -k/m -c/m]; % Fill in 
B = [0; 1/m]; % Fill in
C = [1 0];
p_o=[-2.4, -2.5]; % fill in: Observer pole locations
Ke=place(A', C',p_o); % Ke= transpose of the controller ...
                        % gain matrix of the dual system...
                        % can use place() command 

t_sim=0:.01:10;



%% System calculations
x0=[-2,0]; % released from x=-2 at rest
[t,x]=ode45(@(t,x)smd(x,A,B,u),t_sim,x0);
y=(C*x')'; %([1x2]*[mx2]')' for m time instances
%% Observer calcualtions
x0=[-10,1];
[t_o,x_t] = ode45(@(t_o,x_t)obs(t_o,x_t,A,B,C,u,Ke',y,t),t_sim,x0);

size(x_t)
size(x)
%% Plotting section
% plot the results of x_t and x on the same plots and you can see 
% the performance of the observer. Change the pole locations of the 
% observer and see the difference in behavior of the estimated variables.
% --- Plot results ---
figure(1);
%Displacement
subplot(2,2,1)
plot(t,x(:,1),'b', t,y,'--');
hold on
plot(t,x_t(:,1),'black');
xlabel("time");
ylabel("Displacement (m)");grid on;

%%Velocity
subplot(2,2,2);
plot(t,x(:,2),'b',t,x_t(:,2)); hold on
xlabel("time")
ylabel("Velocity (m/s)");grid on;


%error dynamics
subplot(2,2,3);
plot(t,x(:,1)-x_t(:,1)); hold on
xlabel("time");
ylabel("Error in Displacement estimation")

subplot(2,2,4);
plot(t,x(:,2)-x_t(:,2)); hold on
xlabel("time")
ylabel("Error in velocity estimation(m/s)")

%% System function
% spring mass damper state equations
function dx=smd(x,A,B,u)
% state equations
dx = A*x+B*u;
end

%% observer func

function dx_t = obs(t_o,x_t,A,B,C,u,K_e,y,t)
    y = interp1(t,y,t_o);
    %size(y)
    dx_t = A *x_t + B*u + K_e*(y- C*x_t);
    size(dx_t)
end


