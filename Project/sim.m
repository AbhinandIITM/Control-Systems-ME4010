clc; clear; close all;

%% ========================================================================
% 1. SYSTEM PARAMETERS
% =========================================================================
fprintf('Defining system parameters...\n');
M_w = 0.5;   % Mass of one wheel
M_r = 1.0;   % Mass of robot body
m_b = 0.3;   % Mass of the ball
R = 0.05;    % Radius of the wheel
h = 0.1;     % Height of robot CoM from ball platform
r = 0.02;    % Radius of the ball
g = 9.81;    % Gravity
d = 0.05; 

% Moments of Inertia
I_w = 0.5 * M_w * R^2;        % Inertia of one wheel
I_r = 0.01;                   % Inertia of robot body (given value)
I_b = 0.4 * m_b * r^2;      % Inertia of a solid sphere ball

% Derived parameters
l = h + r; % Height of ball CoM from wheel center

%% ========================================================================
% 2. LINEARIZATION & STATE-SPACE MODEL (A, B)
% =========================================================================
% We linearize the nonlinear equations from `robot_ball_dynamics`
% around the equilibrium point (x=0, theta=0, all velocities=0).
%
% State vector: z = [x; x_dot; theta; theta_dot]
% Input: u = xw_ddot (wheel acceleration)

fprintf('Linearizing nonlinear EOMs to find A and B...\n');

% Linearized Mass Matrix (M_lin * [th_ddot; x_ddot] = F_lin)
% M_lin is constant at equilibrium
M11_lin = I_r + m_b*l^2; % (m_b*x^2 term goes to 0)
M12_lin = m_b*l;
M21_lin = m_b*l;
M22_lin = I_b/r^2 + m_b;

M_lin = [M11_lin, M12_lin;
         M21_lin, M22_lin];
inv_M_lin = inv(M_lin);

% Linearized Forcing Vector (F_lin = F_g*z_pos + F_u*u)
% f1 = (M_r*d*g + g*m_b*l)*theta + (g*m_b)*x - (m_b*l)*u
% f2 = (g*m_b)*theta - (m_b)*u

% F_g matrix (terms multiplying [x; theta])
Fg11 = g*m_b;         % d(f1)/dx
Fg12 = M_r*d*g + g*m_b*l; % d(f1)/dtheta
Fg21 = 0;             % d(f2)/dx
Fg22 = g*m_b;         % d(f2)/dtheta

Fg_mat = [Fg11, Fg12;
          Fg21, Fg22];

% F_u matrix (terms multiplying u)
Fu_vec = [-m_b*l;
          -m_b];
  
% Solve for accelerations in terms of states and input
% [th_ddot; x_ddot] = inv_M_lin * (Fg_mat * [x; theta] + Fu_vec * u)
%
% Let P_mat = inv_M_lin * Fg_mat
% P_mat = [P11, P12]  (row 1: coeffs for th_ddot)
%         [P21, P22]  (row 2: coeffs for x_ddot)
%
% where P11 = coeff of x for th_ddot
%       P12 = coeff of theta for th_ddot
%       P21 = coeff of x for x_ddot
%       P22 = coeff of theta for x_ddot

P_mat = inv_M_lin * Fg_mat;

% Extract coefficients
c_x_x   = P_mat(2,1);  % A(2,1)
c_x_th  = P_mat(2,2);  % A(2,3)
c_th_x  = P_mat(1,1);  % A(4,1)
c_th_th = P_mat(1,2);  % A(4,3)

% G_vec = [b_th; b_x]
G_vec = inv_M_lin * Fu_vec;
b_x   = G_vec(2);      % B(2)
b_th  = G_vec(1);      % B(4)

% Build A and B matrices
% z_dot = A*z + B*u
A = [0, 1,    0,      0;
     c_x_x, 0,  c_x_th, 0;    % x_ddot = c_x_x * x + c_x_th * theta
     0, 0,    0,      1;
     c_th_x, 0, c_th_th, 0];  % th_ddot = c_th_x * x + c_th_th * theta

B = [0;
     b_x;
     0;
     b_th];

fprintf('Newly computed A matrix (with d = %.2f):\n', d); disp(A);
fprintf('Newly computed B matrix (with d = %.2f):\n', d); disp(B);

%% ========================================================================
% 3. LQR CONTROLLER DESIGN
% =========================================================================
fprintf('Designing LQR controller...\n');
% Penalize states: [x, x_dot, theta, theta_dot]
% --- NEW (TUNED FOR DAMPING) ---
Q = diag([70, 60, 90, 50]); % Penalize x_dot and theta_dot
R_i = 0.05; % Penalize control effort

% Check controllability
if rank(ctrb(A, B)) < size(A, 1)
    error('System is NOT controllable. Check parameters or model.');
else
    fprintf('System is controllable.\n');
end

[K, S, E] = lqr(A, B, Q, R_i);
fprintf('New LQR Gain K:\n'); disp(K);

%% ========================================================================
% 4. NONLINEAR SIMULATION
% =========================================================================
fprintf('Running nonlinear simulation with new K and d...\n');
tspan = 0:0.01:30;

% 6-DOF Initial condition: [xw, xw_dot, theta, th_dot, x, x_dot]
x0_6d = [0;             % z(1) = xw
         0;             % z(2) = xw_dot
         deg2rad(6);    % z(3) = theta
         0;             % z(4) = th_dot
         0.01;          % z(5) = x
         0];            % z(6) = x_dot

% Simulate using ode45
% We pass ALL parameters, including the non-zero 'd'
[t, z] = ode45(@(t,z) robot_ball_dynamics(t, z, K, ...
    M_w, M_r, m_b, R, h, r, I_w, I_r, I_b, g, d), tspan, x0_6d);

fprintf('Simulation complete.\n');

%% ========================================================================
% 5. PLOTTING RESULTS (STATIC PLOTS)
% =========================================================================
figure;
set(gcf, 'Color', 'w');

% Plot 1: Ball Position (x)
subplot(2,1,1);
plot(t, z(:,5), 'LineWidth', 1.6);
xlabel('Time [s]');
ylabel('Ball position x [m]');
title('Ball Position (Nonlinear Simulation)');
grid on;

% Plot 2: Robot Tilt Angle (theta)
subplot(2,1,2);
plot(t, rad2deg(z(:,3)), 'LineWidth', 1.6);
xlabel('Time [s]');
ylabel('Tilt angle \theta [deg]');
title('Robot Tilt (Nonlinear Simulation)');
grid on;

%% ========================================================================
% 6. ANIMATION (NEW SECTION)
% =========================================================================
fprintf('Starting animation...\n');
figure;
set(gcf, 'Color', 'w');
ax = gca;
axis equal;
grid on;
hold on; % <-- Make sure hold is on
% Set plot limits (auto-scale x-axis, fixed y-axis)
min_x_anim = min(z(:,1)) - R*2;
max_x_anim = max(z(:,1)) + R*2;
axis([min_x_anim-0.5 max_x_anim+0.5 0 0.5]); 
xlabel('Position (m)');
ylabel('Height (m)');
% --- Define static and dynamic elements ---
% Ground
line(ax, [min_x_anim-5 max_x_anim+5], [0 0], 'Color', 'k', 'LineWidth', 2);
% --- Get initial state (k=1) ---
k=1;
xw_k = z(k,1);
th_k = z(k,3);
x_k = z(k,5);
% --- Calculate initial component positions ---
% Wheel (using rectangle for a circle)
wheel_center = [xw_k, R];
wheel_h = rectangle(ax, 'Position', [wheel_center(1)-R, wheel_center(2)-R, 2*R, 2*R], ...
                    'Curvature', [1 1], 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
% Spoke to show wheel rotation
wheel_angle_k = -xw_k / R; 
spoke_end_x = xw_k + R*sin(wheel_angle_k); % Using sin for x
spoke_end_y = R + R*cos(wheel_angle_k);    % Using cos for y
wheel_spoke_h = line(ax, [xw_k, spoke_end_x], [R, spoke_end_y], 'Color', 'k', 'LineWidth', 1);
% Platform (the line the ball rolls on)
plat_L = 0.3; % Visual length of the platform
Rot_k = [cos(th_k) -sin(th_k); sin(th_k) cos(th_k)];
p1_local = [-plat_L/2; h]; % Platform is at height h=0.1
p2_local = [plat_L/2; h];
p1_world = Rot_k * p1_local + [xw_k; R];
p2_world = Rot_k * p2_local + [xw_k; R];
plat_h = line(ax, [p1_world(1), p2_world(1)], [p1_world(2), p2_world(2)], ...
              'Color', [0.2 0.2 0.8], 'LineWidth', 4); % Thicker platform
              
% Ball
ball_pos_local = [x_k; l]; % l = h+r
ball_pos_world = Rot_k * ball_pos_local + [xw_k; R];
ball_h = rectangle(ax, 'Position', [ball_pos_world(1)-r, ball_pos_world(2)-r, 2*r, 2*r], ...
                   'Curvature', [1 1], 'FaceColor', [1 0.6 0.6], 'EdgeColor', 'r');
title(ax, 'Ball-Bot Animation (Time: 0.00 s)');
% --- Animation Loop ---
% Set target frame rate (e.g., 30 fps)
fps = 30;
sim_step = t(2) - t(1); % Simulation time step (0.01s)
frame_skip = round((1/fps) / sim_step); % Skip frames to match FPS
for k = 1:frame_skip:length(t)
    % Get current state from simulation data
    xw_k = z(k,1);
    th_k = z(k,3);
    x_k = z(k,5);
    
    % Calculate new positions
    Rot_k = [cos(th_k) -sin(th_k); sin(th_k) cos(th_k)];
    
    % 1. Update Wheel
    wheel_center = [xw_k, R];
    set(wheel_h, 'Position', [wheel_center(1)-R, wheel_center(2)-R, 2*R, 2*R]);
    % Update wheel spoke
    wheel_angle_k = -xw_k / R; 
    spoke_end_x = xw_k + R*sin(wheel_angle_k);
    spoke_end_y = R + R*cos(wheel_angle_k);
    set(wheel_spoke_h, 'XData', [xw_k, spoke_end_x], 'YData', [R, spoke_end_y]);
    % 2. Update Platform
    p1_world = Rot_k * p1_local + [xw_k; R];
    p2_world = Rot_k * p2_local + [xw_k; R];
    set(plat_h, 'XData', [p1_world(1), p2_world(1)], 'YData', [p1_world(2), p2_world(2)]);
    
    % 3. Update Ball
    ball_pos_local = [x_k; l]; % l = h+r
    ball_pos_world = Rot_k * ball_pos_local + [xw_k; R];
    set(ball_h, 'Position', [ball_pos_world(1)-r, ball_pos_world(2)-r, 2*r, 2*r]);
    
    % Update plot title with time
    title(ax, sprintf('Ball-Bot Animation (Time: %.2f s)', t(k)));
    
    % Redraw the plot
    drawnow;
end
fprintf('Animation complete.\n');

%% ========================================================================
% NONLINEAR DYNAMICS FUNCTION
% =========================================================================
% This function MUST be in the same file as the main script
function dz = robot_ball_dynamics(t, z, K, ...
    M_w, M_r, m_b, R, h, r, I_w, I_r, I_b, g, d)

    % 1. Unpack the 6-D state vector
    xw      = z(1);
    xw_dot  = z(2);
    theta   = z(3);
    th_dot  = z(4);
    x       = z(5);
    x_dot   = z(6);
    
    % Useful variables
    l = h + r; % Height from ball center to robot CoM
    c_th = cos(theta);
    s_th = sin(theta);

    % 2. Calculate Control Input (LQR)
    % The controller state is [x, x_dot, theta, th_dot]
    state_for_lqr = [x; x_dot; theta; th_dot];
    
    % Control input u is defined as the wheel acceleration
    xw_ddot = -K * state_for_lqr; % This is u
    
    % 3. Build the 2x2 Mass Matrix (Nonlinear)
    % M_mat * [th_ddot; x_ddot] = f_vec
    
    M11 = I_r + m_b*l^2 + m_b*x^2;
    M12 = m_b*l;
    M21 = m_b*l;
    M22 = I_b/r^2 + m_b;
    
    M_mat = [M11, M12; 
             M21, M22];
    
    % 4. Build the Forcing Vector f_vec (Nonlinear)
    % All terms that are NOT th_ddot or x_ddot
    
    % This f1 term correctly includes the non-zero 'd'
    f1 = M_r*d*g*s_th + g*m_b*l*s_th + g*m_b*x*c_th ...
         - m_b*x*s_th*xw_ddot - 2*m_b*x*th_dot*x_dot ...
         - m_b*l*c_th*xw_ddot;
         
    f2 = g*m_b*s_th + m_b*x*th_dot^2 - m_b*c_th*xw_ddot;
    
    f_vec = [f1; f2];

    % 5. Solve for accelerations
    accels = M_mat \ f_vec;
    th_ddot = accels(1);
    x_ddot  = accels(2);

    % 6. Pack the derivatives into dz
    dz = zeros(6,1);
    dz(1) = xw_dot;
    dz(2) = xw_ddot;  % The control input
    dz(3) = th_dot;
    dz(4) = th_ddot;  % Solved
    dz(5) = x_dot;
    dz(6) = x_ddot;   % Solved
end

