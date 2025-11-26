clc; close all; clear;

n1=[1];
d1=[1 3 2];
g=tf(n1,d1)

c=tf([3 8],[1]);
h=tf(5,[1 0]);


% 
% p_g = pole(g);
% z_c = zero(c);
% 
% if all(real(p_g)<0)
%     disp("stable");
% else 
%     disp("unstable");
% end
z = [-0.5]
p = [-1 -2]
k = [2]
sys = zpk(z,p,k)