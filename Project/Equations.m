clc ; clear all 
syms M_w M_r m_b R I_w I_r d h I_b r g T theta(t) x_w(t) x(t) t 
syms x_b y_b 
%%
pos = [cos(theta) sin(theta); -sin(theta) cos(theta)] * [x; h + r] + [x_w; R]; 
% Assign to individual symbols 
x_b = (pos(1)); 
y_b = (pos(2)); 
%%
U_b = simplify(m_b*g*y_b); 
U_r = simplify(M_r*g*(R+d*cos(theta))); 
U_w = simplify(M_w*g*R); 
U = U_b + U_r + U_w  %[output:0c5ba3db]
%%
T_w = simplify(0.5*((M_w*diff(x_w,t)^2)) + 0.5*I_w*(diff(x_w,t)/R)^2); 
T_r = simplify(0.5*(M_r*(diff(x_w/t)^2) + I_r*(diff(theta,t)^2))); 
v_b_x = diff(x_b, t);
v_b_y = diff(y_b, t);

% 2. Calculate translational kinetic energy
T_b_trans = 0.5 * m_b * (v_b_x.^2 + v_b_y.^2);

% 3. Calculate rotational kinetic energy (your term is correct)
T_b_rot = 0.5 * I_b * (diff(x, t)/r)^2;

% 4. Total kinetic energy for the ball
T_b = T_b_trans + T_b_rot;
T_b = simplify(T_b) %[output:3c431c55]
T = T_w + T_r + T_b  %[output:96b97b83]
%%
L = simplify(T - U) %[output:300a5f09]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:0c5ba3db]
%   data: {"dataType":"symbolic","outputData":{"name":"U(t)","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\ng\\,m_b \\,{\\left(x_w \\left(2\\right)+\\sin \\left(\\theta \\left(2\\right)\\right)\\,{\\left(h+r\\right)}+x\\left(2\\right)\\,\\cos \\left(\\theta \\left(2\\right)\\right)\\right)}+M_w \\,R\\,g+\\sigma_1 \\\\\ng\\,m_b \\,{\\left(R+\\cos \\left(\\theta \\left(2\\right)\\right)\\,{\\left(h+r\\right)}-x\\left(2\\right)\\,\\sin \\left(\\theta \\left(2\\right)\\right)\\right)}+M_w \\,R\\,g+\\sigma_1 \n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =M_r \\,g\\,{\\left(R+d\\,\\cos \\left(\\theta \\left(t\\right)\\right)\\right)}\n\\end{array}"}}
%---
%[output:3c431c55]
%   data: {"dataType":"symbolic","outputData":{"name":"T_b(t)","value":"\\left(\\begin{array}{c}\n\\frac{I_b \\,{{\\left(\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)\\right)}}^2 }{2\\,r^2 }\\\\\n\\frac{I_b \\,{{\\left(\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)\\right)}}^2 }{2\\,r^2 }\n\\end{array}\\right)"}}
%---
%[output:96b97b83]
%   data: {"dataType":"symbolic","outputData":{"name":"T(t)","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\n\\sigma_1 \\\\\n\\sigma_1 \n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{I_r \\,{{\\left(\\frac{\\partial }{\\partial t}\\;\\theta \\left(t\\right)\\right)}}^2 }{2}+\\frac{{\\left(M_w \\,R^2 +I_w \\right)}\\,{{\\left(\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)\\right)}}^2 }{2\\,R^2 }+\\frac{I_b \\,{{\\left(\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)\\right)}}^2 }{2\\,r^2 }+\\frac{M_r \\,{{\\left(x_w \\left(t\\right)-t\\,\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)\\right)}}^2 }{t^4 \\,2}\n\\end{array}"}}
%---
%[output:300a5f09]
%   data: {"dataType":"symbolic","outputData":{"name":"L(t)","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\n\\sigma_5 -g\\,m_b \\,{\\left(x_w \\left(2\\right)+\\sin \\left(\\theta \\left(2\\right)\\right)\\,{\\left(h+r\\right)}+x\\left(2\\right)\\,\\cos \\left(\\theta \\left(2\\right)\\right)\\right)}+\\sigma_1 -M_w \\,R\\,g+\\sigma_3 +\\sigma_2 -\\sigma_4 \\\\\n\\sigma_5 -g\\,m_b \\,{\\left(R+\\cos \\left(\\theta \\left(2\\right)\\right)\\,{\\left(h+r\\right)}-x\\left(2\\right)\\,\\sin \\left(\\theta \\left(2\\right)\\right)\\right)}+\\sigma_1 -M_w \\,R\\,g+\\sigma_3 +\\sigma_2 -\\sigma_4 \n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{{\\left(M_w \\,R^2 +I_w \\right)}\\,{{\\left(\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)\\right)}}^2 }{2\\,R^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =\\frac{M_r \\,{{\\left(x_w \\left(t\\right)-t\\,\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)\\right)}}^2 }{t^4 \\,2}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =\\frac{I_b \\,{{\\left(\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)\\right)}}^2 }{2\\,r^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =M_r \\,g\\,{\\left(R+d\\,\\cos \\left(\\theta \\left(t\\right)\\right)\\right)}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_5 =\\frac{I_r \\,{{\\left(\\frac{\\partial }{\\partial t}\\;\\theta \\left(t\\right)\\right)}}^2 }{2}\n\\end{array}"}}
%---
