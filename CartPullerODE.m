m=0.23;
M=2.4;
L=0.36;
g=-9.81;
d=4.4;
%ref=[1;0;pi;0];
t_sim=0:0.001:30;
x0=[0;0;pi/3+0.1;0];
[t,x]=ode45(@(t,x)CtP(x,M,m,L,g,0,d),t_sim,x0);
figure;
hold on;
subplot(2,2,1);
plot(t,rad2deg(x(:,3)),'blue')
ylabel("theta")

subplot(2,2,2);
plot(t,rad2deg(x(:,4)),'blue')
ylabel("theta__dot")

subplot(2,2,3);
plot(t,x(:,1),'blue')
ylabel("x")

subplot(2,2,4);
plot(t,x(:,4),'blue')
ylabel("x__dot")
hold off;

%% CTP
function dx = CtP(x, M, m, l, g, u, b)
    % CtP: cart-pendulum dynamics (point mass at end of rod)
    % states: x(1)=X, x(2)=Xdot, x(3)=theta, x(4)=thetadot
    X      = x(1);
    Xdot   = x(2);
    theta  = x(3);
    thetadot = x(4);

    Den = M + m*(sin(theta))^2;

    % Numerators (derived from solving the 2x2 linear system)
    %EOM not dependant on X.
    Xddot = ( u - b*Xdot + m*l*thetadot^2*sin(theta) - m*g*sin(theta)*cos(theta) ) / Den;

    thetaddot = ( (M + m)*g*sin(theta) + b*Xdot*cos(theta) ...
                  - m*l*thetadot^2*sin(theta)*cos(theta) - u*cos(theta) ) ...
                 / ( l * Den );

    dx = [ Xdot; Xddot; thetadot; thetaddot ];
end
