%% SPRING-MASS-DAMPER SYSTEM WITH STATE OBSERVER
% This script demonstrates state reconstruction using an observer.
% Both noise-free and noisy observer results are plotted together.

clc; clear; close all;

%% ---------------------- System properties ----------------------
m = 10;     % mass (kg)
k = 1000;   % spring constant (N/m)
c = 10;     % damping coefficient (Ns/m)
u = 0;      % external force input (N)

%% ---------------------- System definition ----------------------
A = [0 1; -k/m -c/m];   % state matrix
B = [0; 1/m];           % input matrix
C = [1 0];              % output matrix

% Observer pole locations (must be faster than system poles)
p_o = [-2.4, -2.5];

% Observer gain (dual system approach)
Ke = place(A', C', p_o);

t_sim = 0:0.01:10; % simulation time vector

%% ---------------------- System dynamics ----------------------
x0 = [-2; 0]; % released from x = -2 m at rest
[t, x] = ode45(@(t, x) smd(x, A, B, u), t_sim, x0);
y = (C * x')'; % measured displacement

%% ---------------------- Add measurement noise ----------------------
y_noise = awgn(y, 4); % add 50 dB noise

%% ---------------------- Observer calculations ----------------------
x0_hat = [-10; 1]; % initial estimate

% Noise-free observer
[t_o, x_hat] = ode45(@(t_o, x_t) obs(t_o, x_t, A, B, C, u, Ke', y, t), t_sim, x0_hat);

% Noisy observer
[t_o_noise, x_hat_noise] = ode45(@(t_o, x_t) obs(t_o, x_t, A, B, C, u, Ke', y_noise, t), t_sim, x0_hat);

%% ---------------------- Combined Plots ----------------------
figure('Name', 'Observer Comparison (Noise vs No Noise)', 'NumberTitle', 'off');

% --- Displacement ---
subplot(2,2,1);
plot(t, x(:,1), 'k', 'LineWidth', 1.2); hold on;
plot(t, x_hat(:,1), 'b', 'LineWidth', 1.2);
plot(t, x_hat_noise(:,1), 'r--', 'LineWidth', 1.2);
legend('True displacement', 'Estimated (no noise)', 'Estimated (with noise)');
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Displacement Comparison');
grid on;

% --- Velocity ---
subplot(2,2,2);
plot(t, x(:,2), 'k', 'LineWidth', 1.2); hold on;
plot(t, x_hat(:,2), 'b', 'LineWidth', 1.2);
plot(t, x_hat_noise(:,2), 'r--', 'LineWidth', 1.2);
legend('True velocity', 'Estimated (no noise)', 'Estimated (with noise)');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity Comparison');
grid on;

% --- Error in displacement ---
subplot(2,2,3);
plot(t, x(:,1)-x_hat(:,1), 'b', 'LineWidth', 1.2); hold on;
plot(t, x(:,1)-x_hat_noise(:,1), 'r--', 'LineWidth', 1.2);
legend('Error (no noise)', 'Error (with noise)');
xlabel('Time (s)');
ylabel('Displacement Error (m)');
title('Displacement Estimation Error');
grid on;

% --- Error in velocity ---
subplot(2,2,4);
plot(t, x(:,2)-x_hat(:,2), 'b', 'LineWidth', 1.2); hold on;
plot(t, x(:,2)-x_hat_noise(:,2), 'r--', 'LineWidth', 1.2);
legend('Error (no noise)', 'Error (with noise)');
xlabel('Time (s)');
ylabel('Velocity Error (m/s)');
title('Velocity Estimation Error');
grid on;

sgtitle('Observer Performance: Noise-Free vs Noisy Measurements');

%% ---------------------- Nested Functions ----------------------
% --- System dynamics ---
function dx = smd(x, A, B, u)
    dx = A*x + B*u;
end

% --- Observer dynamics ---
function dx_hat = obs(t_o, x_hat, A, B, C, u, Ke, y, t)
    y_interp = interp1(t, y, t_o);  % interpolate measurement
    dx_hat = A*x_hat + B*u + Ke*(y_interp - C*x_hat);
end
