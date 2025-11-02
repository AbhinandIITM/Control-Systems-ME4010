% System Definition
clear; clc; close all;

s = tf('s');
sys_ol = 16 / (s^3 + 6*s^2 + 5*s + 8);

% Initial Stability Analysis
disp('--- Initial System Analysis ---');
sys_cl = feedback(sys_ol, 1);
poles_cl = pole(sys_cl);
disp('Initial closed-loop poles:');
disp(poles_cl);

figure;
step(sys_cl);
title('Stable Closed-Loop Step Response (Original System)');
grid on;

% Bode and Margin Calculation
disp('--- Calculating Stability Margins ---');
[Gm, Pm, Wcg, Wcp] = margin(sys_ol);

fprintf('Gain Margin (Gm): %.3f (%.2f dB)\n', Gm, 20*log10(Gm));
fprintf('Phase Margin (Pm): %.2f deg\n', Pm);
fprintf('Gain Crossover Freq (Wcg): %.3f rad/s\n', Wcg);
fprintf('Phase Crossover Freq (Wcp): %.3f rad/s\n', Wcp);

figure;
margin(sys_ol);
grid on;
title('Bode Plot with Gain & Phase Margins');

% Time Delay Analysis
pm_deg = Pm;
pm_rad = pm_deg * pi / 180;
time_delay = pm_rad / Wcp; 

s = tf('s');
delay_model = exp(-time_delay * s);
sys_ol_delayed = sys_ol * delay_model;

margin(sys_ol_delayed);

% Critical Gain Analysis
K_crit = Gm; 

sys_cl_crit = feedback(K_crit * sys_ol, 1);
sys_cl_unstable = feedback(1.69 * K_crit * sys_ol, 1); 

% Step Response Plots
figure;
step(sys_cl_crit);
title('Closed-Loop Step Response at Critical Gain (Marginally Stable)');
grid on;

figure;
step(sys_cl_unstable);
title('Closed-Loop Step Response Beyond Critical Gain (Unstable)');
grid on;
