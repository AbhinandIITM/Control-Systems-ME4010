R2 = 1
R3 = 3
nt = 360
theta2 = 2*pi*(0:nt-1)/(nt-1);
% for ii = 1:nt
%     phi(ii) = asin(R2*sin(theta2(ii))/R3)
%     theta3(ii) = pi - phi(ii)
%     R4(ii) = sqrt(R2^2 + R3^2 - 2*R2*R3*cos(pi - theta2(ii) - phi(ii)))
% end
% plot(theta2*180/pi,R4)

phi = asin(R2*sin(theta2)/R3);
theta3 = pi - phi;
R4 = sqrt(R2^2 + R3^2 - 2*R2*R3*cos(pi - theta2 - phi));
plot(theta2*180/pi, R4);