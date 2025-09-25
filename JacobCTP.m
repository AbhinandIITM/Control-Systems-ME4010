% Parameters
M = 1.0;    % cart mass
m = 0.2;    % pendulum mass
l = 0.5;    % pendulum length
g = 9.81;   % gravity
u = 0;      % control force (N)
b = 10;

%%
x0 = [0;0;pi+.1;0] %[output:082e8f85]
%%
A = [0 1 0 0;
    0 -b/M -m*g/M 0;
    0 0 0 1;
    0 -(up*b)/(M*L) 0 
%%
ref = [] %[output:38d15757]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:082e8f85]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0","rows":4,"type":"double","value":[["0"],["0"],["3.2416"],["0"]]}}
%---
%[output:38d15757]
%   data: {"dataType":"error","outputData":{"errorType":"syntax","isTransient":false,"text":"Unsupported use of the '=' operator. To compare values for equality, use '=='. To specify name-value arguments, check that name is a valid identifier with no surrounding quotes."}}
%---
