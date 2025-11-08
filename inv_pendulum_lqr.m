function dx = fcn(x,u)
    dx = zeros(4,1);
    m=1;
    M=5;
    L=2;
    g = -9.81;
    d=1;
    Sx = sin(x(3));
    Cx = cos(x(3));
    D = m*(L^2)*(M+m*(1-Cx^2));

    dx(1,1) = x(2);
    dx(2,1) = (1/D)*(-m^2)
end

function dx = CtP(x,u)
    b = 0;
    
    % CtP: cart-pendulum dynamics (point mass at end of rod)
    % states: x(1)=X, x(2)=Xdot, x(3)=theta, x(4)=thetadot
    X      = x(1);
    Xdot   = x(2);
    theta  = x(3);
    thetadot = x(4);

    %Den = M + m*(sin(theta))^2;
    Den = M + m*(sin(theta))^2;
    % Numerators (derived from solving the 2x2 linear system)
    %EOM not dependant on X.
    % Xddot = ( u - b*Xdot + m*l*thetadot^2*sin(theta) - m*g*sin(theta)*cos(theta) ) / Den;
    % 
    % thetaddot = ( (M + m)*g*sin(theta) + b*Xdot*cos(theta) ...
    %               - m*l*thetadot^2*sin(theta)*cos(theta) - u*cos(theta) ) ...
    %              / ( l * Den );
    % Xddot = ( u - b*Xdot + m*g*sin(theta)*cos(theta) + m*l*(thetadot^2)*sin(theta) ) / Den;
    % 
    % thetaddot = ( (M + m)*g*sin(theta) + cos(theta)*( u - b*Xdot + m*l*(thetadot^2)*sin(theta) ) ) ...
    %         / ( l * Den );
    Xddot = ( u - b*Xdot + m*l*thetadot^2*sin(theta) + m*g*sin(theta)*cos(theta) ) / Den;
     
    thetaddot = ( (M + m)*g*sin(theta) - b*Xdot*cos(theta) ...
                   + m*l*thetadot^2*sin(theta)*cos(theta) + u*cos(theta) ) ...
                  / ( l * Den );

    dx = [ Xdot; Xddot; thetadot; thetaddot ];
end

%[appendix]{"version":"1.0"}
%---
