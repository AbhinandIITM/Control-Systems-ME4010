x0 = [10,1];
m = 10 ; k = 1000; u = 0; c = 1;
t_sim = 0:0.01:10;
[t,x] = ode45(@(t,x) smd(x,m,k,c,u), t_sim,x0);
plot(t, x(:,1),'LineWidth',1.5); xlabel('Time (s)'); ylabel('displacement'); grid on
%% smd
function dx = smd(x, m, k, c, u)
    dx(1,1) = x(2);
    dx(2,1) = -k/m * x(1) - c/m * x(2) + u/m;
    size(dx)
end


