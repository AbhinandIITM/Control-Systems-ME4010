clc ; clear all 
syms M_w M_r m_b R I_w I_r d h I_b r g T theta(t) x_w(t) x(t) t 
syms x_b y_b 
%%
pos = [cos(theta) sin(theta); -sin(theta) cos(theta)] * [x; h + r] + [x_w; R]; 
% Assign to individual symbols 
x_b = simplify(pos(1)); 
y_b = simplify(pos(2)); 
%%
U_b = simplify(m_b*g*y_b); 
U_r = simplify(M_r*g*(R+d*cos(theta))); 
U_w = simplify(M_w*g*R); 
U = U_b + U_r + U_w  %[output:0c5ba3db]
%%
T_w = simplify(0.5*((M_w*diff(x_w,t)^2)) + 0.5*I_w*(diff(x_w,t)/R)^2); 
T_r = simplify(0.5*(M_r*(diff(x_w/t)^2) + I_r*(diff(theta,t)^2))); 
T_b = simplify(0.5*(m_b*((diff(x,t)*cos(theta) - x*sin(theta)*diff(theta,t) + (R+r)*cos(theta) + diff(x_w,t))^2+(-diff(x,t)*sin(theta) -x*cos(theta)*diff(theta,t) - (h+r)*sin(theta))^2) + I_b * (diff(x,t)/r)^2)); 
T = T_w + T_r + T_b  %[output:15b31016]
%%
L = simplify(T - U) %[output:4ef2aa14]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:0c5ba3db]
%   data: {"dataType":"symbolic","outputData":{"name":"U(t)","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\ng\\,m_b \\,{\\left(x_w \\left(2\\right)+\\sin \\left(\\theta \\left(2\\right)\\right)\\,{\\left(h+r\\right)}+x\\left(2\\right)\\,\\cos \\left(\\theta \\left(2\\right)\\right)\\right)}+M_w \\,R\\,g+\\sigma_1 \\\\\ng\\,m_b \\,{\\left(R+\\cos \\left(\\theta \\left(2\\right)\\right)\\,{\\left(h+r\\right)}-x\\left(2\\right)\\,\\sin \\left(\\theta \\left(2\\right)\\right)\\right)}+M_w \\,R\\,g+\\sigma_1 \n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =M_r \\,g\\,{\\left(R+d\\,\\cos \\left(\\theta \\left(t\\right)\\right)\\right)}\n\\end{array}"}}
%---
%[output:15b31016]
%   data: {"dataType":"symbolic","outputData":{"name":"T(t)","value":"\\frac{I_r \\,{{\\left(\\frac{\\partial }{\\partial t}\\;\\theta \\left(t\\right)\\right)}}^2 }{2}+\\frac{m_b \\,{\\left({{\\left(\\cos \\left(\\theta \\left(t\\right)\\right)\\,\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)+\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)+\\cos \\left(\\theta \\left(t\\right)\\right)\\,{\\left(R+r\\right)}-\\sin \\left(\\theta \\left(t\\right)\\right)\\,x\\left(t\\right)\\,\\frac{\\partial }{\\partial t}\\;\\theta \\left(t\\right)\\right)}}^2 +{{\\left(\\sin \\left(\\theta \\left(t\\right)\\right)\\,\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)+\\sin \\left(\\theta \\left(t\\right)\\right)\\,{\\left(h+r\\right)}+\\cos \\left(\\theta \\left(t\\right)\\right)\\,x\\left(t\\right)\\,\\frac{\\partial }{\\partial t}\\;\\theta \\left(t\\right)\\right)}}^2 \\right)}}{2}+\\frac{{\\left(M_w \\,R^2 +I_w \\right)}\\,{{\\left(\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)\\right)}}^2 }{2\\,R^2 }+\\frac{I_b \\,{{\\left(\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)\\right)}}^2 }{2\\,r^2 }+\\frac{M_r \\,{{\\left(x_w \\left(t\\right)-t\\,\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)\\right)}}^2 }{t^4 \\,2}"}}
%---
%[output:4ef2aa14]
%   data: {"dataType":"symbolic","outputData":{"name":"L(t)","value":"\\begin{array}{l}\n\\left(\\begin{array}{c}\n\\sigma_6 +\\sigma_1 -g\\,m_b \\,{\\left(x_w \\left(2\\right)+\\sin \\left(\\theta \\left(2\\right)\\right)\\,{\\left(h+r\\right)}+x\\left(2\\right)\\,\\cos \\left(\\theta \\left(2\\right)\\right)\\right)}+\\sigma_2 -M_w \\,R\\,g+\\sigma_4 +\\sigma_3 -\\sigma_5 \\\\\n\\sigma_6 +\\sigma_1 -g\\,m_b \\,{\\left(R+\\cos \\left(\\theta \\left(2\\right)\\right)\\,{\\left(h+r\\right)}-x\\left(2\\right)\\,\\sin \\left(\\theta \\left(2\\right)\\right)\\right)}+\\sigma_2 -M_w \\,R\\,g+\\sigma_4 +\\sigma_3 -\\sigma_5 \n\\end{array}\\right)\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{m_b \\,{\\left({{\\left(\\cos \\left(\\theta \\left(t\\right)\\right)\\,\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)+\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)+\\cos \\left(\\theta \\left(t\\right)\\right)\\,{\\left(R+r\\right)}-\\sin \\left(\\theta \\left(t\\right)\\right)\\,x\\left(t\\right)\\,\\frac{\\partial }{\\partial t}\\;\\theta \\left(t\\right)\\right)}}^2 +{{\\left(\\sin \\left(\\theta \\left(t\\right)\\right)\\,\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)+\\sin \\left(\\theta \\left(t\\right)\\right)\\,{\\left(h+r\\right)}+\\cos \\left(\\theta \\left(t\\right)\\right)\\,x\\left(t\\right)\\,\\frac{\\partial }{\\partial t}\\;\\theta \\left(t\\right)\\right)}}^2 \\right)}}{2}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_2 =\\frac{{\\left(M_w \\,R^2 +I_w \\right)}\\,{{\\left(\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)\\right)}}^2 }{2\\,R^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_3 =\\frac{M_r \\,{{\\left(x_w \\left(t\\right)-t\\,\\frac{\\partial }{\\partial t}\\;x_w \\left(t\\right)\\right)}}^2 }{t^4 \\,2}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_4 =\\frac{I_b \\,{{\\left(\\frac{\\partial }{\\partial t}\\;x\\left(t\\right)\\right)}}^2 }{2\\,r^2 }\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_5 =M_r \\,g\\,{\\left(R+d\\,\\cos \\left(\\theta \\left(t\\right)\\right)\\right)}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_6 =\\frac{I_r \\,{{\\left(\\frac{\\partial }{\\partial t}\\;\\theta \\left(t\\right)\\right)}}^2 }{2}\n\\end{array}"}}
%---
