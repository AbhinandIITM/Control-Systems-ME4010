% Put this in a new file (or run in the command window)
% Parameters (replace with your real values)
I_w = 0.01;    % wheel inertia
I_r = 0.02;    % rotor/pivot inertia
I_b = 0.03;    % body translational inertia (or equivalent)
R   = 0.05;    % wheel radius
r   = 0.01;    % small radius offset (as in your eqns)
M_r = 0.5;     % rotor mass (or equivalent)
M_w = 0.6;     % wheel mass
m_b = 1.2;     % body mass
h   = 0.2;     % geometry constant
d   = 0.1;     % (if used)
g   = 9.81;

params.I_w = I_w; params.I_r = I_r; params.I_b = I_b;
params.R = R; params.r = r; params.M_r = M_r; params.M_w = M_w;
params.m_b = m_b; params.h = h; params.g = g; params.d = d;

% time span and initial condition (states: [x_w; x_wdot; theta; thetadot; x; xdot])
tspan = 0:0.001:50;
% choose initial x_w = 0, x_wdot = 0, theta = pi/3, thetadot = 0, x = 0, xdot = 0
x0 = [0; 0; 0+deg2rad(5); 0; 0.1; 0];

[t, X] = ode45(@(t,X) threeDOF_dynamics(t,X,params), tspan, x0); %[output:68240a29]

% plot theta and x (same layout as your cartpuller script)
figure; %[output:6fee1f84]
subplot(2,2,1); plot(t, rad2deg(X(:,3))); ylabel('theta (deg)'); %[output:6fee1f84]
subplot(2,2,2); plot(t, rad2deg(X(:,4))); ylabel('theta dot (deg/s)'); %[output:6fee1f84]
subplot(2,2,3); plot(t, X(:,5)); ylabel('x (m)'); %[output:6fee1f84]
subplot(2,2,4); plot(t, X(:,6)); ylabel('x dot (m/s)'); %[output:6fee1f84]

% ---------------- ODE function (save as threeDOF_dynamics.m) ----------------
function dX = threeDOF_dynamics(~, X, p)
% State ordering:
% X(1)= x_w, X(2)= x_wdot, X(3)= theta, X(4)= thetadot, X(5)= x, X(6)= xdot

x_w    = X(1); x_wdot = X(2);
theta  = X(3); thdot  = X(4);
x      = X(5); xdot   = X(6);

% parameters
I_w = p.I_w; I_r = p.I_r; I_b = p.I_b;
R   = p.R;   r   = p.r;
M_r = p.M_r; M_w = p.M_w;
m_b = p.m_b; h   = p.h; g = p.g;

% --- Build mass/inertia matrix M (3x3) so that M * [xdd_w; thdd; xdd] + f = 0
% Note: entries follow the coefficient of the corresponding acceleration in the
% three scalar EOMs you provided. Check signs carefully with your printed derivation.

% Row 1 (wheel equation): coefficients multiplying [xdd_w, thdd, xdd]
M11 = I_w + R*(M_r + M_w) - R;  % I_w* xdd_w  + R(M_r+M_w)*xdd_w  - R*xdd_w (from -R[... + xdd_w])
M12 = - m_b*(R + r)*sin(theta)*x;   % coefficient of thdd : -m_b(R+r) sinθ * x
M13 = - R * x * cos(theta);         % coefficient of xdd : -R * (x cosθ)

% Row 2 (theta equation): coefficients multiplying [xdd_w, thdd, xdd]
% I_r*thdd  + (R m_b sin(2θ)/2)*(xdd + r*xdd_w)  - h*m_b*sin^2(θ)*xdd  + m_b*(r-x)*sinθ * xdd_w  ...
M21 = (R*m_b*sin(2*theta)/2)*r + m_b*(r - x)*sin(theta); % combined coefficient multiplying xdd_w
M22 = I_r;                                              % coefficient of thdd
M23 = (R*m_b*sin(2*theta)/2) - h*m_b*(sin(theta))^2;    % coefficient of xdd

% Row 3 (body x equation): coefficients multiplying [xdd_w, thdd, xdd]
M31 = m_b*(cos(theta) - x);          % coefficient of xdd_w : ( + m_b cosθ  - m_b x )
M32 = - R*m_b*sin(2*theta);          % coefficient of thdd
M33 = I_b + m_b;                      % coefficient of xdd  (I_b + m_b*1)

M = [ M11, M12, M13;
      M21, M22, M23;
      M31, M32, M33 ];

% --- Build f vector (all remaining terms NOT multiplied by accelerations)
% These include Coriolis-like terms, centrifugal, gravity, etc. Move them to RHS as f.
% Form: M*acc + f = 0  =>  acc = -M\f

% Row 1 residuals (from your eqn1):
% f1 collects terms: -m_b(R+r) sinθ*(2 xdot thdot + x thdot^2 cotθ) - R*( x*thdot^2 )
f1 = - ( - m_b*(R + r)*sin(theta)*( 2*xdot*thdot + x*thdot^2*cot(theta) ) ...
         - R*( x*(thdot)^2 ) );
% NOTE: signs: original eqn had "- m_b(...) - R[ ... ] = 0". After moving acceleration terms to M,
% the leftover terms enter f (see above). If simulation shows sign-flips, flip sign of f1/terms.

% Row 2 residuals (theta eqn)
% collect everything not multiplied by thdd or xdd or xdd_w:
% - M_r*d*g*sin(theta)  + R*m_b*sin(2θ)* xdot + R*m_b*r*sin(2θ)*thdot^2
% - g*h*m_b*sin(theta)  + m_b*( 2*x*thdot*xdot ) 
f2 = - ( - M_r * p.d * g * sin(theta) ...
         + R*m_b*sin(2*theta)* xdot ...
         + R*m_b*r*sin(2*theta)*thdot^2 ...
         - g*h*m_b*sin(theta) ...
         + m_b*( 2*x*thdot*xdot ) );

% Row 3 residuals (x eqn)
% terms: - h*m_b*sin(2θ)*thdot^2 - g*m_b*sinθ + ( - m_b*x* xdd_w term is in M31 ) + Coriolis ...
f3 = - ( - h*m_b*sin(2*theta)*thdot^2 - g*m_b*sin(theta) ...
         + ( - m_b*x* sin(theta) * 0 ) ... % placeholder if you had extra mixed terms
         + ( 2*m_b*x*thdot*xdot ) );

% NOTE: The above f1,f2,f3 are assembled from the text you gave. Because the algebra in the
% text is dense, you should double-check each sign/term against your symbolic derivation.
% If you want, paste your original LaTeX/symbolic EOMs and I'll translate them term-by-term.

f = [ f1; f2; f3 ];

% Solve for accelerations
acc = - M \ f;   % acc = [ xdd_w; thdd; xdd ]

xdd_w = acc(1);
thdd  = acc(2);
xdd   = acc(3);

% assemble state derivative
dX = zeros(6,1);
dX(1) = x_wdot;
dX(2) = xdd_w;
dX(3) = thdot;
dX(4) = thdd;
dX(5) = xdot;
dX(6) = xdd;
end


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:68240a29]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Failure at t=7.180165e-01.  Unable to meet integration tolerances without reducing the step size below the smallest value allowed (1.776357e-15) at time t."}}
%---
%[output:6fee1f84]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAgoAAAE6CAYAAAB3U3plAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnX+IXtd554\/DsmTKZnFn7RKEa0\/cHW391zRxDdNJwJjAKt1lhjaL0YwKFcP8Me7GkMVWNTPZgBDU0owamQrFxfPH7KAsq5FwQ9YVbfBql9SQqDIpJPPHkiC1ykhdDyZxFVFCFJalszzXOW\/OXN33fe+59557z4\/PC8YzmvPjeT7POed+3\/PrPrS3t7en+EAAAhCAAAQgAIECAg8hFGgXEIAABCAAAQj0I4BQoG1AAAIQgAAEINCXAEKBxgEBCEAAAhCAAEKBNgABCEAAAhCAgD0BZhTsmZEDAhCAAAQgkAwBhEIyocZRCEAAAhCAgD0BhII9M3JAAAIQgAAEkiGAUEgm1DgKAQjYELh\/\/7565ZVX1NGjR9X4+LhNVtJCICoCCIWowokzEIBAEYGbN2+qU6dOqbNnz6rR0dFekrW1NbW+vp79fvHiRTU5Odn72+XLl9X58+fV5uYmQoFmlTQBhELS4cd5CMRPQETC\/Py8evTRR9XGxkZPKFy\/fl2JUJB\/u3HjRu9nERKS58KFC+ojH\/mI+uxnP4tQiL+Z4OEAAggFmgcEIBAtAZkVWFlZUYuLi0qEgSkURCTIZ2lpSckyw\/Lyspqbm1MTExO9JYe33npLHTp0CKEQbQvBsTIEEAplKJEGAhDwkoA8\/G\/fvq0OHz6c2Se\/b21tqdXVVTUyMqK2t7fVwYMHs\/\/r2QOZMdDCYGpqKstr\/v6JT3wim4HY3d3Nypyenu6V5yUEjIKAYwIIBceAKR4CEHBHQD\/g79y5oz796U9nSwZaJJi1mssMplCQGQS9L0GExNjYWE90SP4vf\/nLzCi4Cx8lB0IAoRBIoDATAhDoT8BcRihKhVCg9UCgOgGEQnV25IQABDomUHdGoWjpQS9jdOwa1UPAGwIIBW9CgSEQgIAtgWF7FHR5+RkF+XdzqcHczGgekbS1h\/QQiJEAQiHGqOITBCCwj0CRUBh0PBJ8EIDALwggFGgNEIBA9ASKhIKeVeh34VL0UHAQAiUJIBRKgiIZBCDQHgF9\/4HUePr06X0nEdqzgpogAAEhgFCgHUAAAl4RkCOOx48fV2fOnMns0j\/zvgWvwoQxCRFIWijwrSWhlo6rwRCQfnnt2rXefQhF9xsE4wyGQiACAskKBb61RNB6cSFKAvk7EYbdkRAlBJyCgEcEkhUKfGvxqBV2YMraWz\/Iar1z92dWtT8++uEs\/dKhj1nlI3F5AvkZBOmrOzs72TsZzM\/ntr5XvlBSJkVA99N+TtN\/7ZpDskKh7LeWR154o0f0Qz\/9h0K6H\/rp+0aaD37+8Pf\/3C4SpLYm8LNfn1H\/9EuPZPl+8f9\/1fu5X4H741Uc00HG\/NMv9a9Dly1t5Z+9\/\/1W28GtW7esGfqYoYxQePLJJ9VPPnXcR\/OxyQMC0kerfPQ4YptX+v2\/\/B\/7hawuI4Z+mbRQMO917\/et5cC\/+0\/q85\/\/fBbzft8+8\/\/+9z\/+4Ftq2W+rWv3+6i9\/8G1Vfpf\/6qpeGUxDaqTyLV8zk\/+\/88476sDBiUKO5jcG4aZ\/f23uKds+3nh67ccbX39b\/b9H\/k2vfLFRbL3yuY83XqcUGFq8+0EoI+JD9zVk+7H9wZYrfX7trR1199Xn9v0xZFamIwiFn791rp9QaCPQ+Wlw\/bAUwTFMbOQfmFpomP+vKzhsnmrmw94US2XEkymYfHrw2\/g\/KK2w+ebf3lPf+rt7PUHYhHCQNtrGtxdz82\/ezyaPMOb7YtFmxjb6ZVNxp5x0CEy\/9p3MWf1lIJZ2mrRQkIDqdc9+G6ZCCHQ\/oSH+lXlAN9WNB4kWqcOHb\/tN+dpEOTpu8k1EC7vvfvG3Khftoq3evXtXLSwsZK9pHiQGtIg4cOCA2tzcVHWOMpbZaOzC18rgyQgBg8DoS99Qf\/4fP64+9a8fjmaWL1mhUOZbi8SeAYkxoC0CevoyE7CHxqyXnppuqyISTp48qU6cOKHk1cxlPjrPuXPnyiTvm2bY0eWmfa1lLJkhYBD4jT\/6a\/XJX3s4+2IUSztNViiU+daCUKD\/d0VAdvRvffu9bO9F2VmGWAalMsxT8rUMD9L4Q8DcrxBLO01WKEizGvatBaHgT+dL1RI96JQRDK4HJb0MIct1Bw8e7C1JTExMqI2NjdKzDk3E0rWvTdhIGekSkOUH2dgYSztNWiiUacaxBLqMr6Txl0C\/XdWmxa7bqrmPR0T2pUuXMoFw9erVwnsOXNJ07atL2yk7fgIIhfhjvM9DBqTEAu6xu3JiYuZPv9N3\/4LLtmrOJsgMwvLyspKNizK70O\/NjC5RuvTVpd2UnQYBhEIace55yYCUWMADcFcPQnlTXbbVomWH2dnZ7K2OCIUAGg0mtkpANjTK0Xa5iOn9159vtW4XlbH0MISqy8HXRUApM34CstFR7mLIb3J02Vbv37+fzSJMTU2pJ554Qh07dqx3DLKLdzG49DX+FoSHrgnozcgIBdekPSmfAcmTQGDGPgIyqyDHr+ae+Whrs19yUmh+fl7t7u6qxcXFbNlBRIL8vrq6qkZGRlqLEv2yNdRUVJGA9FGEQkV4oWVjQAotYmnYq1+IZF5ilVJbTcnXNFp0nF7Ku4JYeogztvu8YkBKIMiBupjfq9B0W9X7Eh5\/\/PHWZwyGhaRpX4fVx98hUIUAQqEKtQDzMCAFGLRETM4vP7hqq+aSg6Bt8r0OVUPlyteq9pAPAkUEEAqJtAsGpEQCHaCb+eWHttqqeVFZE+92qIK+LV+r2EYeCGgCCIVE2gIDUiKBDtDN\/CVMXbRVfRrizp07rd7O2IWvATYRTO6YAEKh4wC0VT0DUlukqacKAXOfQkptNSVfq7QL8vhBAKHgRxycW8GA5BwxFdQg0JZQMF83PchcfWyyhkulstIvS2EiUccEEAodB6Ct6hmQ2iJNPVUImBsaXbdV8\/0O+rXTWkDILY0zMzP7rnau4k\/ZPK59LWsH6SAwiABCoWb7MDdE5YvyYVe1tokBqWagye6UgFwVK5cuLR36mNM31ZlXOE9OTu7zybzC+caNG9klTF\/72tec+k2\/dIqXwhsigFCoANKcvhwkBrSI6GpHtekaA1KFQJOlNQLTr31HySuo5eIll20VodBaSKkoIgIIBctgykBz8uRJdeLEidLvrdd5zp07Z1lbc8ldDr7NWUlJqRKQI5Ly8pkrn\/u4U6EgfIctPcgLoiTNtWvXlOs+S79MtcWH5TdCIax4VbaWAakyOjK2QECOSG59+73sBVFttFVZZjhy5Mg+zy5evKhkOUJEwvnz53svi3Lpfhu+urSfstMggFBII86tDL6JoMRNBwREJMiswt1Xn0uqrSIUHDQmimycAEKhBlLfjloNcoUBqUagydoKAX1EMqW2mpKvrTQiKnFCAKFQE+uw9c42j1ohFGoGk+ydEmhTKJinlWTJ4fbt29meBF4z3WkToHJPCSAUagTGtx3UCIUawSRr5wTaEgpy7HF3d1cdP35cvfjii2ppaUlNTEy0dneCCZoZhc6bHQaUIJB\/w2uJLF4meWhvb2+vbcsQCm0Tp76YCbQhFMw+e\/DgQbWwsJAJBdnEaN6joC9ics0boeCaMOU3QQChUJPisKWHNo9aMaNQM5hk75QAQqFT\/FQOgb4EEAoNNA5fjlohFBoIJkV0RqANoSDO6TsSzKUHPbsgVziLuG\/rw4xCW6Sppw4BhEIdegHlZUAKKFiJmirXOHd5j0IXV67TLxNt7IG5jVAILGBVzWVAqkqOfG0R0O97WP+DT6tbt261VW2n9dAvO8VP5SUJIBRKghqUzJejViw9NBBMiuiMgH7fw9f\/879HKHQWBSqGwIMEEAo1W4VPR60QCjWDSfZOCbgSCmUvRhPn5ZjkxsZG6fe41AXGjEJdguRvgwBCoQZl345aIRRqBJOsnRPQL4b632f\/g9MZhUEnlfRRybZgIBTaIk09dQggFGrQQyjUgEdWCOQItCEUyt590sQ9Cvfv388ucbpy5UrP08XFxezeBv1BKNANQiCAUKgZJZ+OWjGjUDOYZO+UgH6D5D\/+l99zNqPQplCQul5++WX1hS98QY2PjxeyRSh02uSovCQBhEJJUIOSFd2j0MVRK4RCA8GkiM4I6DdIPvzfF5wJBXGuzCVpTUC4efOmOnXqlDp79mzfPQ8IhSZIU4ZrAggF14Q9KZ8ByZNAYMZAAjIguRYKYsCgS9KaClGZK6Hpl03RphyXBBAKLul6VDYDkkfBwJS+BNoSCm2EwDw2LfUVnaigX7YRCeqoSwChYEnQ56NWLD1YBpPk3hGISSjoo9P61dX53wW+CAX5pHLBlHcNDoNKEUAolMLUP5FPR60QCjWDSfbOCbgQCiLuT548qU6cOFH6fgSd59y5c6WYyH6E+fn57PXV8rl48WL2RkrzI2nk\/RJnzpzpbW5kRqEUXhJ1TAChUCMAbe6grmFmlpUBqS5B8rdBwIVQELvNmcCih7j2TS8XHDhwQG1ubvY9rVCFRdHmRvplFZLkaZsAQqEGcYRCDXhkhUABAVdCwaxKlgDW19cL+Td1WknfoTA1NZW9jVL\/LgKEexRo+qERQCjUjFhbR61qmsmMQl2A5G+FgLwY6v3\/+afq\/\/zVf22lPpeV5C9cmp6eVnq\/gq6XGQWXEaDspgggFBog2cZRq7pmMiDVJUj+NgiIUPjhtctq9y\/\/pI3qOq+Dftl5CDCgBAGEQglIMSRhQIohivH7gFCIP8Z4GB4BhEJ4MatkMUKhEjYytUwAodAycKqDQAkCCIUSkMwkbRy1sjSpVHKEQilMJOqYAEKh4wBQPQQKCCAUKjQLF0etzJ3YRUezzL\/nj3eZN8D127WNUKgQaLK0TkCEwnvf\/V\/qvT\/7Yut1d1Eh\/bIL6tRpSwChYEssl76Jo1b6DZR6R3T+JIV5Z\/yNGzeU1LmxsZFdHmNe4iKm5S900eYyINUMNNlbITD92nfUt7\/xF86Egm9HmumXrTQrKqlJAKFQE6CL7Pkb3EQYyEfOX+sjV3Nzc9nNb3mRIWnHxsays9vmhwHJRaQos2kCCIWmiVIeBOoTQCjUZ9h4CaZQeOyxx9Ty8rLKX9yifzdFhBiS\/50ZhcbDQ4EOCbgSCvkXNPVzYXFxcd+FSA5dzYpGwLsmTPlNEEAoNEGx4TLMl8dI0SIU9AyCFgN61iA\/gyAD4s7OzgODHQNSw0GiuMYJSBv9yaeOqw\/99P1Olh4ad6hEgfTLEpBI0jkBhELnIdhvgDzoz58\/37tnPr\/UgFDwLGCY0ygBVzMKjRrZYGEIhQZhUpQzAggFZ2iLCx70lrm8SJAS+t0Zz9JDy4GjulYIfG7re+qrX\/2qsxkF7UTRUkRT73mwAYVQsKFF2q4IIBS6Ip+rt+idETqJubxQtJnRXGpgM6MnAcWMSgREKLzx9bfVD7\/yQqX8ZTL59H4WhEKZiJGmawIIhZoRyM8QmMVNTEz0jjEOqkaOPx47dqzva205HlkzSGQPhoBrocDxyGCaAoZ6RAChUCMY5rLAzMxMb9PhwYMH1cLCQrahUI4wDvv0u4vBvFiJC5eGUeTvMRBAKMQQRXyIjQBCoUZE899OzGl\/mQXY2tp64LWyNaqrlZUpzlr4yNwSAddCQdxg6aGlYFJNNAREKLw295Sae+ajQfv00N7e3l7bHuSFgnk00VwukBsUu\/4gFLqOAPWXIdCGUNBiYWVlZZ9JbGYsEyHSpEgAoVAz6uYFR6Y4uHr1qrp27RozCjX5kj0tAm0JBV+oIuB9iQR2DCKAUKjZPvLHF\/VegqIXO9WsqlZ2BqRa+MjcEgGEQkugqQYCFgQQChawQk6KUAg5eunY7kIomG97HUay7EmlYeWU\/Tv9siwp0nVJAKFQg75vR60GucKAVCPQZG2NgAuhkDd+0GbGsieVmgJCv2yKJOW4JIBQqEEXoVADHlkhUEDAtVDwrc8iFOgGIRBAKFSIkq9vomNGoUIwyeIVAYSCV+HAGAhkBBAKNRrCoG8nNYp1kpVvLk6wUmjDBFwLBTGXexQaDhrFRU8AoRB9iD9wEKGQSKADd7MNoSCI5CjzkSNH9tEyb0JtCyP9si3S1FOHAEKhDr2f5zWXImSwuX37tld3KCAUGggyRbRCoC2h0IozJSpBKJSARJLOCSAUaoZA7k3Y3d1Vx48fVy+++GL2fgc5YrW8vKzkLgX53YcPA5IPUcCGYQQQCsMI8XcItE8AoVCDublHIf8iKK5wrgGWrMkSQCgkG3oc95gAQqFGcBAKNeCRFQIFBNbe+oH646\/9jXr\/9eeT4MNMXxJhDt5JhELNEMr+BHmng7n0oGcXZmdn1eHDh2vW0Ex2BqRmOFKKWwIIBbd8KR0CVQggFKpQy+Up2kHdxZvoBrmCUGgg0BThnABCwTliKoCANQGEgjWyMDMgFMKMW2pWtyEUbt68qebn57NNyPkP73pIrcXhbxkCCIUylCJIg1CIIIgJuOBaKJhve52ZmclOJ83Nzan8ZuS2UNMv2yJNPXUIiFBYOjSmlg59rE4xned9aG9vb69zKzw2gAHJ4+BgWo+Aa6GQv01VjjePjY1le4lkCXFra0utrq6qkZERq6hoASKiY3JyspfXvGOlaDmSfmmFmcQdEUAo1ATv0zTmIFcYkGoGmuytEGhbKMiDfGdnJ7vvpOqRZi0Srly5oszbHWVskE3OZ86cydjpn8fHx3ss6ZetNCsqqUkAoVADoB4gfLpYqZ87DEg1Ak3W1gi4FgriiMwiyCcvDq5evWp9o6r+ovD000+rO3fuZGXqGQV9IkrPUJizFxoo\/bK1pkVFNQggFGrA46VQNeCRFQIFBNoQCuY+BVlykAf4+vp6dpPq5uamMr\/xDwvSu+++myWRpYqFhYV9QsEUJHmBglAYRpa\/+0QAoVAjGv3WJWsU6Swr31ycoaXgBgm0IRQaNLdXVNGXhvwMgrnMYQoF+fnWrVsuzKJMCDRCAKFQE2PVdc2a1VpnRyhYIyNDBwRcC4VBs4B1+nIdoYBI6KChUaUVAYSCFS6l9ICwvb09NGfbZ7IHGYRQGBouEnhAwGehoJcoBNP09PS+0xH9hIKk1S+Gyy9FyN\/olx40OkwYSgChMBRRHAkYkOKIY+xeuBIK5jHFQQwXFxcrvfG1SCjklxrYzBh7643XP4RCjdi6msasYVLfrAgFF1Qp0wWBR154w9lLoVxtQC4ql+ORLloHZXZBAKFQgzpCoQY8skKgDwEZlO6++lxQfPqNBVy4FFQYMbYPARHvf\/i7v8nNjDYtxPU0po0tZdMyo1CWFOm6JtCGUDD7sFySdPv2bes7FJrgRL9sgiJluCaAUKhB2NU0Zg2TWHpwAY8yWyXgWijIXgF5IZT5anjZeCzvfWj78jSEQqtNi8oqEkAoVAQXWjYGpNAilq69LoWCKe7zL4KqczyyarTol1XJka9NAgiFNml3WBcDUofwqdqKAELBCheJIeCcAELBOWI\/KkAo+BEHrBhOwKVQkNr1OxjMpQc9uzA7O5u9SbKtD\/2yLdLUU4cAQqEOvYDyMiAFFKzETXUtFASvLDMcOXJkH+mi10C7DgX90jVhym+CAEKhCYoBlMGAFECQMDEjkFJbTclXmne4BBAK4cbOynIGJCtcJO6QQEptNSVfO2xSVF2TAEKhJsBQsjMghRIp7HTdVuXGxPn5+eyIZP7T9vtZXPtKa4JAEwQQCk1QDKAMBqQAgoSJzpce9Kvh274voV9o6Zc0+hAIIBRCiFIDNjIgNQCRIloh4LKt+nZJmktfWwkWlSRBAKGQRJjT2iCWSEijddPlw1PPKMzNzanJycnOGbr0tXPnMCAaAgiFaEI52BEGpEQCHYGbrttqFzcwsvQQQcNM2AWEQiLBdz34JoIRN1sg0HRb1csN29vbQ61nM+NQRCRIkMCv\/P7r6vnffla9NvdU0N4\/tLe3txe0Bz83Xk+NTk1N7bshTl5ks76+nqWSt92Z06bDXmUreZoefGNgjQ9+Ekiprabkq5+tDavKEEAolKHUYhr90DdviTOnSm\/cuKFENGxsbKjR0VElR73kKtozZ85kVuqfx8fH91nNgNRiEKmqFoGU2mpKvtZqFGTulABCoVP8+yuXKdKXX35Z3bt3T5l3zoswkM\/S0pLKb8bS99avrq6qkZGRTESMjY09cF89A5JHgcaUgQRcttVBpx662Lvg0leaGQSaIoBQaIpkA+XIQ\/6ZZ55Rb775ptJLD\/mliPzvpogQE\/K\/a7MYkBoIEEW0QsBlW0UotBJCKomMAELBk4DKEsKFCxfUSy+9pE6ePPmAUDCPc5mzBvkZBJlh2NnZyWYfzI8MvnwgEAqBW7duNWqquY9nUMGLi4sP9J1GDckVRr90SZeymyLwk08dZzNjUzCrliOzBK+88oo6evSoeuyxx9Ty8nLjQqGqbeSDQEwEfLtwKSa2+AIB3wkEc+ohf8+8nGCQz9tvv71vD0LTSw++BxD7IAABCEAAAi4JBCMUiiCYRx\/Nv+tpUHN5oWgzo7nU0G8zo0v4lA0BCEAAAhDwnUDQQsGEW3SPQhPHI30PIPZBAAIQgAAEXBKIWigIuLoXLrmET9kQgAAEIAAB3wlEIxSaBj1IYDRdV53yzL0b09PTSt8LMWypRl4XvLm5qfIXTNWxpUpeG\/t1+f1u4axSf508Nrbnr0PO3xJax47U8obQN23ahumPD\/3Sxnbf+qTYY2M\/\/bLc6IFQKOA0aMmiHNZ2UpkPzJmZmX2nPvIW5C+Ykt8vXbrUu6myHYv312Jjv5mz6BbOtu23sV2nlYeAHL81bwXtWqi1za1ufSH0TZu24Vu\/tLHdtz4p9tjYT78s3xsRCgWsBt3oWB6t+5T5B44MoltbWwNnFbRVPjysqtjf7xZO97T312Bju6Q9deqUOnv2bHZ9OJ\/qBELomzZtI0+i635ZxXZf+qSeTTCv4x80JtIvy\/dDhEKO1bAbHcujdZ8yf3WuzVW6XQ9IQqeK\/UW3cLon\/WANNrbnvzV2YW8MdYbSN23ahm9CoYrtvvRJ2zGFfll+VEAo9BEK\/W50LI\/Wfcq8WrZRyNK5d3d3S80+uPLE1v5+t3C6sm9QuTa261s\/pbx+bzLtwofQ6swfcRb7fTzWbNM28jHoul\/a2u5Tn9RCwZxVHTQm0i\/LjwAIhQSFgnSQ8+fPd76Z0WZQGnQLZ\/nm3lxKG9v1ngq9gVHyHjt2rHP+zdFop6TYhYIP\/dKmXfvWJ6sIhZWVFUW\/HN5\/EQqJLT34MBhp5DbTnJK23y2cw5t58ylsbM9PcfpyaqN5Km5LjHnpwZd+adOufeuTdZce6Jf9+y9CoYDNoBsd3Q6FdqXnp9WGbWb04aSD6aGN\/cNu4bQjVz+1je35uBR9M65vURolhNA3bdqGRM2nfmlju299Ulja2E+\/LD9mIBQKWIVwBEvMtjkK5ON0t439Zph8UP42tudfqGSz6bR8V04jZQh906Zt+NYvbWz3rU\/ajon0y\/JjBkKhD6sQLnXRCnp+fj7bmJi\/cMl8dXY\/9d\/1xT+DLkfp9+pvH4SCDXtJa17s4sOlOuWHCP9ShtA3y7ZrH\/tlWdt9FAr0Szf9FaHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIggFCIIow4AQEIQAACEHBDAKHghiulQgACEIAABKIgkKxQuHv3rlpYWFDb29tZIKenp9Xq6qoaGRmJIrA4AQHXBG7evKnm5+fV7u5uVtXi4qJaWlrqVWv+vah\/ra2tqfX19Sz9xYsX1eTkZC\/v5cuX1crKSvb76dOn1eHDh127Q\/kQgEAfAkkKhfv376vl5WU1NTWVDUD69wMHDuwb6Gg1EIBAMQEttEUYyANe\/z47O7uvT0kfm5mZ2dffpMTr168rEQobGxvqxo0bvZ9HR0eVCIzjx4+rM2fOZJXrn8fHxwkHBCDQAYEkhUIRZ\/kGc+3aNWYVOmiEVBkHAXnwy0fEg\/mwlwe8CIOtra1e\/zLTaqE+NzeXiY58X5S0Y2NjzCrE0UzwIkACCIWfBw2hEGDrxWSvCJgPf3PGQGYJzN9lea9oRk\/P8JnliIP5371yGmMgkAABhIJSD0ybJhB3XIRAowT0foQvfelL2axAfgZB\/n7q1Cl19uzZbB+QCAU9g6DFgJ41yM8giIjf2dlhWbDRiFEYBMoTSF4o6GlPQVa0mfHJJ58sT5OUEOiQwK1btzqpXe9PEIGgNzO6Fgr0y05CTaWWBLrqk5ZmDk2etFAYJhKEngxIoQY7ZNtDZ9+F\/aMvfUPdffW5oZ2+yQRFIkHKd730QNtuMop2ZYXMvm3b267PLpLlUycrFMqedIgl0OWbBClDJdC2UMifdDC5mUsNeo9CfjOjXmoo2sxoLjUUbWakX4baStOyO5Z2mqxQkMFHzn8PuzshlkCn1T3T9LZNoTBMaJtHkF0cj6RfptnGQ\/M6lnaapFDIX7akG9\/ExER2rlu+AelPLIEOrYNhrz2BNoVC\/rIlba15sZLLC5fol\/btgxztE4ilnSYpFGyaSyyBtvGZtGESaFModE2Iftl1BKi\/DIFY2ilCYUi0Ywl0mUZNmrAJIBTCjh\/Wx0cglucHQgGhEF\/vTNQjhEKigcdtbwkgFLwNTbOGxRIY6xKkAAAVEklEQVToZqlQmo8EEAo+RgWbUiYQy\/ODGQVmFFLux1H5jlCIKpw4EwEBhEIEQSzjQiyBLuMracImgFAIO35YHx+BWJ4fzCgwoxBf70zUI4RCooHHbW8JIBS8DU2zhsUS6GapUJqPBBAKPkYFm1ImEMvzgxkFZhRS7sdR+Y5QiCqcOBMBAYRCBEEs40IsgS7jK2nCJoBQCDt+WB8fgVieH8woMKMQX+9M1COEQqKBx21vCSAUvA1Ns4bFEuhmqVCajwQQCj5GBZtSJhDL84MZBWYUUu7HUfmOUIgqnDgTAQGEQgRBLONCLIEu4ytpwiaAUAg7flgfH4FYnh\/MKDCjEF\/vTNSjqkLh8uXLamVlpZDa6dOn1eHDh70jGssA7B1YDGqUQCztFKGAUGi0Y1BYdwRshMLdu3fVwsKC2t7eVoPEgBYRBw4cUJubm2p8fLw7B42aYxmAvYCJEc4I\/Mrvv66e\/+1n1WtzTzmro42CEQoIhTbaGXW0QKCsUBCRcPLkSXXixAk1OjpayjKd59y5c6XSu06EUHBNmPKbIIBQaIJiAGUwIAUQJEzMCJQVCjHgol\/GEMX4fUAoxB\/jzEMGpEQCHYGbCIUIgogLURFAKEQVzv7OIBQSCXQEblYVCnq\/wtLSkjp48GBv78LExITa2NgovTzRJkL6ZZu0qasqAYRCVXKB5WNACixgCZtbVSisra1l1EQoyObFS5cuZQLh6tWramdnJ\/t33z70S98igj1FBBAKibQLBqREAh2Bm1WEgjmbIDMIy8vLSk44iDi4fv26EhHh46wC\/TKCBpuACyIUXv69Q2rp0MeC9pZTD0PCx4AUdPtOyvi6QkEvO8zOzmZ3JyAUkmo+OOuAwCMvvKH+8Hd\/E6HggK1XRSIUvAoHxgwgUEUo3L9\/P5tFmJqaUk888YQ6duxY774Ec0nCN\/D0S98igj1FBBAKibQLBqREAh2Bm1WEgrh98+ZNNT8\/r3Z3d9Xi4mK27CAiQX5fXV1VIyMj3tGhX3oXEgwqIIBQSKRZMCAlEugI3KwqFEJ0nX4ZYtTSsxmhkEjMGZASCXQEbpYVCnoD4+OPP+7tjMGwcNAvhxHi7z4QQCj4EIUWbGBAagEyVTRCoKxQ0JWZSw7yb76+AKoIDv2ykSZDIY4JIBQcA\/aleAYkXyKBHcMI2AqFfHnmWyR9ewlU3lb65bDWwN99IIBQ8CEKLdjAgNQCZKpohIAMSu+\/\/nwjZenTEHfu3OEehUaIUkiKBBAKiUQdoZBIoCNws6pQyC9BmCh8vcaZfhlBg03ABYRCAkEWFxmQEgl0BG5WEQp65kDfxhgKBvplKJFK206EgkfxN9dW82bV3aDFgORRoDFlIIEqQsG8wnlycjIYwvTLYEKVtKEIhY7Drwe47e3tgbu1tYioujmLAanjQFN9KQJrb\/1A\/fHX\/sZ6j4KeUZibm1MIhVKoSQSB0gQQCqVRNZ9QRMLJkyfViRMnSr8CV+c5d+6clUEIBStcJO6IQFWhIOb6\/E6Hfjjplx01NKq1IoBQsMIVbmIGpHBjl5LldYQCmxlTain42iYBhEKbtDusC6HQIXyqLk2gqlAwXwolb4wM5UO\/DCVSaduJUPAo\/i6\/ETEgeRRoTOlLoKpQYDMjjQoC7gggFNyxtSrZ9fEuhIJVOEjcEYGqQoHNjB0FjGqTIIBQ8CTMrr8RIRQ8CTRmDCRQVShIoTIjd+rUKXX27NnSm4O7Dgf9susIUH8ZAgiFMpRaSOP6GxEDUgtBpIraBKoKBfOYcZER3MxYOzQUkDABhIJHwXd5vAuh4FGgMaUvgapCIVSk9MtQI5eW3QgFj+LNZkaPgoEpnRBAKHSCnUohMJAAQsGTBtLE8a61tTU1Njamio6H8c3Fk0BjxkACNkKh6QvL+i3\/mQJ+enpara6uqpGRkZ4f0u\/W19ez3y9evLjvZkjzWvaia9jpl3SIEAggFDyJUt3NjHqw6vdOCAYkTwKNGY0JBSnI3JuQf0ibFQ27Al2LhCtXrux72JsCfmZmRi0vL6upqameGDeXC2\/cuKGkH25sbGSbKUVgHD9+XJ05cyYzRf88Pj7eM41+SYcIgQBCwZMoVd3MqAfKhx9+OPPkM5\/5DDMKnsQUM+wJ2Mwo5Es3v9nn\/zbopWp6xuDpp59Wd+7cUUtLS71ZAfNhLw94EQZbW1u9WQWpUz6SJ9+HRZxcu3ZtX9r8jB9Cwb6NkKN9AgiF9pn3rbHK8S4RCj\/+8Y+VvCwq\/23HrIgByaNAY0pfAnWEQlWs7777bpZVlhMWFhb2CYX8BmPzd0lv9rn88qEpIqT8\/O\/yb\/TLqlEjX5sEEApt0h5QV93jXcP2ODAgeRJozBhIoAuhoA0qWv7LzyCYYl4LBfONleY+ofyeIZlh2NnZyYSI\/tAv6RAhEBh96Rvqtbmn1NwzHw3B3L42PrS3t7cXtAc1jUco1ARI9s4JyEPzZ78+o\/7v45+0fs10E8Z3JRTE9lu3bjXhAmVAoHEC0i\/v\/c4GQqFxsiULbHLXNkKhJHSSeU3AxxkFc4MiSw9eNx+Mc0SAGQVHYMsW2\/SubXNHtmkDU5xlI0K6Lgn4JhTy+4aKNjPqDYpFmxnNpYai48v0yy5bG3WXJYBQKEuqhXRVd22LacwotBAgqnBOoKpQGHS8uOyNp0VlcDzSecipIAACCIUAglTGRIRCGUqk8Z2ACIWz\/+0t9cOvvGBlqiuhIEZw4ZJVKEgcIQGEQoRBLXKJKc5EAh24m5\/b+p564+tvlxYK5s2Hg1xfXFzcd9rAF0z0S18igR2DCCAUEmkfDEiJBDpwN22Fgna37s2mXWGjX3ZFnnptCCAUbGgFnJYBKeDgJWR6VaEQKiL6ZaiRS8tuhEIi8WZASiTQgbtZVygULUUMur65a1z0y64jQP1lCCAUylBqIY1sRpT\/5GUyLj4MSC6oUmbTBOoIBREJly5d6r2USWzTSxKzs7OF70Bp2n7b8uiXtsRI3wUBhEIX1AvqHDagXbhwQckrbqsKCQYkTwKNGQMJVBUKTZx66CI09MsuqFOnLQGEgi0xh+nlvPeRI0cyQaDfea+PZj366KP7vinZmsGAZEuM9F0QQCh0QZ06ITCYAELBsxaivxn96Ec\/UiIOtre3VRNHuxAKngUacwoJVBUKUhhLDzQqCLghgFBww7VWqeYFL01txEIo1AoJmVsiUEcoaLGwsrKyz9qm+pALBPRLF1Qps2kCCIWmidYsT1\/jLLMIzz777ANLEVWLZ0CqSo58bRKoKxTatLWJuuiXTVCkDNcEEAquCZcs31xy2NzcVOPj41lO\/e\/y88bGBpsZS\/IkWZgEEAphxg2r4yaAUPAkviIIrly5oo4ePVpoEacePAkUZjglYCMUzDevDjNqYmKiltAeVn7VvzOjUJUc+dokIELh7qvPtVmlk7oe2tvb23NSciSFMiBFEsjI3bARCnkUgzYzLi0tqcnJSe\/o0S+9CwkGFRBAKCTSLBiQEgl04G5WFQrcoxB44DHfawIIBa\/D05xxCIXmWFKSOwIIBXdsKRkCVQkgFKqSCywfQiGwgCVqblWhILi4RyHRRoPbzgkgFJwj9qMChIIfccCKwQTqCAUpWd9uatZy8eJFL\/cniI30S3pECAQQCiFEqQEbGZAagEgRzgnUFQrODWy4Avplw0ApzgkBhIITrP4VyoDkX0yw6EECCAVaBQT8I4BQ8C8mTixCKDjBSqENE0AoNAyU4iDQAAGEQgMQQygCoRBClLARoUAbgIB\/BBAK\/sXEiUUIBSdYKbRhAgiFhoFSHAQaIIBQaABiCEUgFEKIEjZWFQr3799X8t\/o6GhQEOmXQYUrWWMRComEngEpkUAH7ub0a99R77zzjvrhV16w8kTfzDg7O6sOHz78QN6670qxMsYiMf3SAhZJOyOAUOgMfbsVMyC1y5vaqhEQofDtb\/yFeu\/PvmhdgL5DYXp6Wq2urqqRkRF18+ZNNT8\/rx599FFeCmVNlAwQ+IAAQiGRloBQSCTQAbu59tYP1Df\/9l5loSCum69rF3Gwvb2tFhcXlbwUyscP\/dLHqGBTngBCIZE2wYCUSKADdlMGI\/n88zvfqjSjoF3Xswi7u7vq9OnThUsRvmCiX\/oSCewYRAChkEj7YEBKJNCBuilLDvL51t\/dUx\/+\/ptq9y\/\/pJIna2tran19PZtFePbZZ9WRI0eUuRRRqVCHmeiXDuFSdGMEEAqNofS7IAYkv+OTunUyEL0295Sae+ajld5\/YC45bG5uqvHx8Qyp\/nf5eWNjw7tTEfTL1Ft+GP4jFMKIU20rGZBqI6QAhwR+44\/+Wn33i7+V1VClrYoguHLlijp69GihlZx6cBg8io6eAEIh+hB\/4GCVwTcRNLjZIYGtb7+n5O6EpUNjaunQx5Jrq\/TLDhsfVZcmgFAojSrshAxIYccvRuvllMPaWzvq7qvP7XMvpbaakq8xtuFUfEIoJBJpBqREAh2AmzKLICJBPnq5wTQ7pbaakq8BNE1M7EMAoZBI02BASiTQHrs5TCBo01Nqqyn56nHTxLQhBBAKiTQRBqREAu2hmzJ7ICJBPp\/8tYez0w2DPim11ZR89bBpYlJJAgiFkqBCT8aAFHoEw7JfNijeufuz7F6Ex0c\/XLjE0M+jlNpqSr6G1YKx1iSAUEikPTAgJRLoDtzUouDvf\/yzTBzIx1YcmGan1FZT8rWDpkmVDRFAKDQE0vdiGJB8j5Df9okYkI+eJTCtFVFQZkmhrIcptdWUfC0bf9L5RwCh4F9MnFjEgOQEa9CFmg9\/ccScEcgLAfldxIB8hu0xqAslpbaakq912wX5uyOAUOiOfas1MyC1irv1ymTDoJ721\/\/v9+AX42QWQD6\/+ssf7v3sWgCUhZJSW03J17LxJ51\/BBAK\/sXEiUUMSE6wOim06kPf1we\/LaSU2mpKvtq2A9L7QwCh4E8snFrCgOQU777CzQe9\/MH8hm\/+XmSR+U3f\/OYv\/66vOG7Pk25qSqmtpuRrN62JWpsggFBogmIAZTAg\/SJI+lZA\/QDXf8n\/LlP3wx7s+dDrB73+dp\/qw75Ol0iprabka502Qd5uCSAUuuXfWu0S6NQ\/5kNcs5A1ev3J\/z2lb\/E+tY2YHp6XL19WKysrGd7Tp0+rw4cP70Mdk68+tSFsaZYAQqFZnp2UNmwwEqMYkDoJDZVWIBBLW71586Y6fvy4OnPmTEZB\/zw+Pt6jEouvFcJMloAIIBQCClaRqWUGI4RC4EFOzPxYHp4i4K9du6ZWV1fVyMiIWltbU2NjY\/tmFWLxNbEmmpy7CIXAQ15mMEIoBB7kxMyP5eEpwkA+S0tL2f\/zv9MvE2vYAbuLUAg4eEWDT9FgxIAUeJATMz8moWDOIIio39nZ6QkH+mViDTsgd\/XJLXlXi7wK\/pEX3lDvv\/58QB4Um\/rQ3t7eXvBeVHAgP51ZNBjpAalC8WSBQCcEbt261Um9TVZapm+KKLr3OxtNVktZEKhF4EM\/fT\/L\/6Gf\/oP6F988o\/7x334wM4ZQqIW128xlBqNuLaR2CKRJoMzSQ5pk8DokAjK7EMsdLknPKEijG7QOGlKjxFYIxEIgP7tXtJkxFl\/xAwIhEEhWKDAYhdA8sTFFAmVPJKXIBp8h0AWBZIUCg1EXzY06IVCOQJk7TsqVRCoIQKAugWSFgoBjMKrbfMgPAQhAAAKxE0haKAwKrqyLrq+vZ0kuXryoJicnvWwLMjMyPz+vdnd31fT0dO+SmiJjTZ8OHDigNjc3lXnbXRcO2tiv7bt\/\/75aXl5WU1NTD1zt26YPNrbfvXtXLSwsqO3tbe\/bVJsMq9QVQt+0aRu+9Usb233rk2KPjf30y3I9EKFQwOn69evZJS8bGxvqxo0bvZ9HR0fLUW0plfnAnJmZGfjwzF8wJb9funQp87Erv2zsN5HqmaCidwC0hF7Z2K7TijiTzbPmslfXQq0tXk3VE0LftGkbvvVLG9t965Nij4399MvyvRKhUMDKPJ6lG9Pc3Jx3swr5B44MoltbWwNnFbS7Pjysqtgv3wBefvllde\/ePTU7O9vZjIKN7ZL21KlT6uzZs52JsvJDgt8pQ+ibNm0jT7vrflnFdl\/6pJ5NMN8NMmhMpF+W7+sIhRyr\/LS2L9PcRSE1v13JrED+90HNoOsBSWyrYr88KJ555hn15ptvdrr0YGN7\/ltj+e5JSpNAKH3Tpm34JhSq2O5Ln7QdU+iX5ccXhEIfoWDOIPh6jjuvlm0Usvgk+xr0i3fKN5nmUtraL\/5duHBBvfTSS+rkyZOdCwVz9mYQe30UV8iFsO+luQg3W1LR7J6PfdO2XZuUuu6Xtrb71Ce1UKBfNtvvpDSEQoJCQR5c58+f73wzo82gJA+JV155RR09elQ99thjnW9mtLFd76nQm2Il77Fjxzrn3\/xw4rbE2IWCD\/3Spl371ierCIWVlZXeZnX6Zf\/+i1BIbOnBh8FII7eZ5pS0b7\/9drYZ0IflIBvb81OcPtjv9pHupvSYlx586Zc27dq3Pll36YF+iVCwGrnM6UzfNzOam+SGbWb04aSDGYj8dP0g+80jZGYZi4uL+94qaBXoGoltbM\/75XObqoGklawh9E2btiHQfOqXNrb71ieFpY399MvyXZYZhQJWIRzBErNtjgL5OK1mY78ZJh+Uv43t+qy2zIbIfRw2m07Ld+U0UobQN23ahm\/90sZ23\/qk7ZhIvyw\/ZiAU+rAK4VIXraD7Xbhkvs+in\/rv+jKpQZej9Hv1tw9CwYa9pDUvdvHlsqvyw4RfKUPom2XbtY\/9sqztPgoF+qWbvopQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IYBQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IYBQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IYBQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IYBQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IYBQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IYBQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IYBQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IYBQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IYBQcMOVUiEAAQhAAAJREEAoRBFGnIAABCAAAQi4IfD\/AfvY\/OdurHC4AAAAAElFTkSuQmCC","height":251,"width":417}}
%---
