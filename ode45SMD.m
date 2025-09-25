x0 = [10,0];
m = 10 ; k = 1000; u = 0;
t_sim = 0:0.01:10;
[t,x] = ode45(@(t,x) smd(x,m,k,c,u), t_sim,x0);
plot(t, x(:,1),'LineWidth',1.5); xlabel('Time (s)'); ylabel('displacement'); grid on