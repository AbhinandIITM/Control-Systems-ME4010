% clear; clc;
%%

m1 = 3; c1 = 2; k1 = 4;
m2 = 5; c2 = 5; k2 = 12;
u = 0;

%%

A = [0               1               0               0;
    -(k1+k2)/m1   -(c1+c2)/m1       k2/m1          c2/m1;
     0               0               0               1;
     k2/m2          c2/m2          -k2/m2         -c2/m2];

B = [0; 0; 0; 1/m2];
C = [1 0 0 0];
%%
eig(A) %[output:56baf6ad]
%%

poles_o = [-6 -6.5 -7 -8];

K = place(A', C', poles_o)';   
%%
t_sim = 0:0.01:10;
x0 = [1; 0; -1; 0];
%%
[t, x_true] = ode45(@(t, x) smd(x, A, B, u), t_sim, x0);

%%
y_true = (C * x_true')';

y_noise = awgn(y_true, 10, 'measured');
%%
xhat0 = [0; 0; 0; 0];
%%

[t_o, x_hat] = ode45(@(t_o, xhat) ...
                     fso_dyn(t_o, xhat, A, B, K, u, y_true, t), ...
                     t_sim, xhat0);

[t_o2, x_hat_noise] = ode45(@(t_o, xhat) ...
                     fso_dyn(t_o, xhat, A, B, K, u, y_noise, t), ...
                     t_sim, xhat0);

%%

figure('Name','Full-State Observer','NumberTitle','off');

for i = 1:4
    subplot(2,2,i);
    plot(t, x_true(:,i), 'k','LineWidth',1.2); hold on;
    plot(t_o, x_hat(:,i), 'b','LineWidth',1.2);
    plot(t_o2, x_hat_noise(:,i), 'r--','LineWidth',1.2);
    xlabel('Time (s)');
    ylabel(['State x',num2str(i)]);
    legend('True','Estimated','Estimated (noise)');
    grid on;
    title(['State x',num2str(i),' Estimation']);
end

sgtitle('Full-State Observer');

%% -------- Functions --------
function dx = smd(x, A, B, u)
    dx = A*x + B*u;
end

function dxhat = fso_dyn(t, xhat, A, B, L, u, y_vec, t_vec)
    y_interp = interp1(t_vec, y_vec, t);   % measurement at this time
    dxhat = A*xhat + B*u + L*(y_interp - [1 0 0 0]*xhat);
end
%%


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":30.9}
%---
%[output:56baf6ad]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"ans","rows":4,"type":"complex","value":[["-1.5594 + 2.2046i"],["-1.5594 - 2.2046i"],["-0.1073 + 0.6537i"],["-0.1073 - 0.6537i"]]}}
%---
