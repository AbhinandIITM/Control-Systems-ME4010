function dx = CtP(x,M,m,l,g,u,b)
    % States
    X      = x(1);      % cart position
    Xdot   = x(2);      % cart velocity
    theta  = x(3);      % pendulum angle
    thetadot = x(4);    % angular velocity
    
    
    denom = m*(l^2)*(M + m*sin(theta)^2);
    
    Xddot = (u*m*(l^2) + m*l*(thetadot^2)*sin(theta) -b*Xdot -(m^3)*(l^2)*g*cos(theta)*sin(theta))/ denom;
    thetaddot = ((m+M)*m*g*l*sin(theta) - m*l*cos(theta)*(m*l*(thetadot^2)*sin(thetadot) - b*theta) -u*l*cos(theta))/denom;
    dx = [Xdot; Xddot; thetadot; thetaddot];
end