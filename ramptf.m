% sys = tf([1],[1 1 1])
% t = 0:0.01:10
% ramp = (t.^2)/2
% [y t] = lsim(sys,ramp,t)
% plot(t,y)

syms s t

% Transfer function: 1 / (s^2 + s + 1)
G = 1/(s^2 + s + 1);
U = 1/s^3;
Y = G * U;
y_t = ilaplace(Y, s, t);
simplify(y_t)
fplot(y_t, [0 10])
grid on
xlabel('Time (s)')
ylabel('y(t)')
title('Symbolic response to parabolic input')

