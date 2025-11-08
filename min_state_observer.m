
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
eig(A)
%%
x0 = [1; 0; -1; 0];
xhat0 = [0; 0; 0];
%%
A_aa = [0]  %[output:4ab6c872]
A_ab = [1 0 0]  %[output:8d8de9fc]
A_ba = [-(k2+k1)/m1 0 k2/m2]'  %[output:3de400c2]
A_bb = [-(c1+c2)/m1 k2/m1 c2/m1; 0 0 1; c2/m2 -k2/m2 -c2/m2]  %[output:76084649]
B_a = [0]  %[output:0c2cf9ba]
B_b = [0; 0; 1/m2]  %[output:36c7d452]
C_a = [1]  %[output:8994f6e4]
C_b = [0 0 0] %[output:463d8005]
%%
min_order_poles = [-6 -6.5 -7 ];
K_e = place(A_bb',A_ab',min_order_poles)' %[output:82178fe2]
%%
F_hat = B_b - K_e*B_a %[output:184931c2]
A_hat = A_bb - K_e*A_ab %[output:2a47b13b]
B_hat = A_ba - K_e*A_aa + A_hat*K_e %[output:1b43e133]
C_hat = [0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1] %[output:37d7b671]
D_hat = [1 ; K_e] %[output:596fa4d3]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:4ab6c872]
%   data: {"dataType":"textualVariable","outputData":{"name":"A_aa","value":"0"}}
%---
%[output:8d8de9fc]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"A_ab","rows":1,"type":"double","value":[["1","0","0"]]}}
%---
%[output:3de400c2]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"A_ba","rows":3,"type":"double","value":[["-5.3333"],["0"],["2.4000"]]}}
%---
%[output:76084649]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"A_bb","rows":3,"type":"double","value":[["-2.3333","4.0000","1.6667"],["0","0","1.0000"],["1.0000","-2.4000","-1.0000"]]}}
%---
%[output:0c2cf9ba]
%   data: {"dataType":"textualVariable","outputData":{"name":"B_a","value":"0"}}
%---
%[output:36c7d452]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"B_b","rows":3,"type":"double","value":[["0"],["0"],["0.2000"]]}}
%---
%[output:8994f6e4]
%   data: {"dataType":"textualVariable","outputData":{"name":"C_a","value":"1"}}
%---
%[output:463d8005]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"C_b","rows":1,"type":"double","value":[["0","0","0"]]}}
%---
%[output:82178fe2]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"K_e","rows":3,"type":"double","value":[["16.1667"],["2.5875"],["58.1500"]]}}
%---
%[output:184931c2]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"F_hat","rows":3,"type":"double","value":[["0"],["0"],["0.2000"]]}}
%---
%[output:2a47b13b]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"A_hat","rows":3,"type":"double","value":[["-18.5000","4.0000","1.6667"],["-2.5875","0","1.0000"],["-57.1500","-2.4000","-1.0000"]]}}
%---
%[output:1b43e133]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"B_hat","rows":3,"type":"double","value":[["-197.1500"],["16.3187"],["-985.8850"]]}}
%---
%[output:37d7b671]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"C_hat","rows":4,"type":"double","value":[["0","0","0"],["1","0","0"],["0","1","0"],["0","0","1"]]}}
%---
%[output:596fa4d3]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"D_hat","rows":4,"type":"double","value":[["1.0000"],["16.1667"],["2.5875"],["58.1500"]]}}
%---
