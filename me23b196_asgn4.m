%%
% RLOCUS of initial transfer function and finding dominant poles and
% corresponding gain value using click feature
G = zpk([-5],[-2 -4 -7 -9],1)
rlocus(G)
sgrid(0.7,[])
%Dominant Poles are -2.81+2.87i and gain value K=88.9 from the plot
%%
% Q1 Finding the compensator pole for given compensator zero
% Given Parameters

zeta = 0.7; % Desired damping ratio

Ts = 1; % Desired settling time

desired_comp_zero = -4.5; % The compensator zero you want to test

% Original System's Poles and Zeros

plant_poles = [-2, -4, -7, -9];

plant_zeros = [-5];

sys=zpk(plant_zeros,plant_poles,88.9);

step(sys)

% 1. Calculate the Desired Dominant Pole Location

sigma = 4 / Ts;

omega_n = sigma / zeta;

omega_d = omega_n * sqrt(1 - zeta^2);

s_desired = -sigma + 1i * omega_d;

disp(' ');

disp(['Desired dominant pole location: ', num2str(s_desired)]);

% 2. Calculate the Angle Deficit (using the Angle Criterion)

% Angle from plant poles and zeros

angle_poles_sum = sum(angle(s_desired - plant_poles));

angle_zeros_sum = sum(angle(s_desired - plant_zeros));

% Angle from the new compensator zero

angle_comp_zero = angle(s_desired - desired_comp_zero);

% The sum of angles from all poles and zeros must equal -180 degrees (or 180).

% We can use -180 for convenience.

angle_total_sum_rad = angle_zeros_sum + angle_comp_zero - angle_poles_sum;

angle_total_sum_deg = rad2deg(angle_total_sum_rad);

% The required angle from the compensator pole to satisfy the angle criterion

% is the remaining angle to reach -180 degrees.

angle_comp_pole_required_deg = 180 + angle_total_sum_deg;

angle_comp_pole_required_rad = deg2rad(angle_comp_pole_required_deg);

disp(['Required angle from compensator pole: ', num2str(angle_comp_pole_required_deg), ' degrees']);

% 3. Calculate the Compensator Pole Location from the Required Angle

% The angle is defined by atan(imag(s_desired) / (real(s_desired) - real(p_c)))

% We can rearrange this to solve for the real part of the compensator pole.

comp_pole_real_part = real(s_desired) - imag(s_desired) / tan(angle_comp_pole_required_rad);

% The compensator pole is a real pole on the real axis

compensator_pole = comp_pole_real_part;

disp(['Calculated compensator pole location: ', num2str(compensator_pole)])
%%
%step response of uncompensated system for given k value
G = zpk([-5],[-2 -4 -7 -9],[88.9]);
sys=G;
% desired damping ratio
rlocus(G)
sgrid(0.7, [])  % damping ratio line
step(sys);
    title('Unit Step Response For Uncompensated System');
    xlabel('Time (seconds)');
    ylabel('Output');
   grid on;
%%
%Testing for different compensator zero locations and corresponding step responses
zeros=[-5 desired_comp_zero]
poles=[-2 -4 -7 -9 compensator_pole]

% 2. Define the Gain Calculation Function
% (This function should be saved in a separate file named 'calculate_gain_K.m')
function K = calculate_gain_K(zeros, poles, s_desired)
% CALCULATE_GAIN_K Calculates the gain K required for the root locus to pass through s_desired.
% K = (Product of pole magnitudes) / (Product of zero magnitudes) at s_desired

    % Calculate the product of the magnitudes of the distances from s_desired to ALL poles
    pole_magnitudes = abs(s_desired - poles);
    product_of_pole_magnitudes = prod(pole_magnitudes);

    % Calculate the product of the magnitudes of the distances from s_desired to ALL zeros
    zero_magnitudes = abs(s_desired - zeros);
    product_of_zero_magnitudes = prod(zero_magnitudes);

    % Apply the magnitude criterion: K = Poles / Zeros
    K = product_of_pole_magnitudes / product_of_zero_magnitudes;
end

% 3. Execute the Function and Display Result
K_calculated = calculate_gain_K(zeros, poles, s_desired);
%Step response of compensated system
G1 = zpk([-5 desired_comp_zero],[-2 -4 -7 -9 compensator_pole],[K_calculated]);
%This gain value is calculated manually using K=magnitude of product of
%poles/magnitude of product of zeros formula at new dominant pole location
sys1=G1;
rlocus(G1)
step(sys1)
  title('Unit Step Response for original compensated system');
    xlabel('Time (seconds)');
    ylabel('Output');
    grid on;
    
  
%%
% Q2
% Given Parameters
zeta = 0.738;          % Desired damping ratio
Ts = 1;              % Desired settling time
desired_comp_pole = -8.82; % The compensator pole you are using

% Original System's Poles and Zeros
plant_poles = [-2, -4, -7, -9 -8.82];
plant_zeros = [-5 -4.5];

% 1. Calculate the Desired Dominant Pole Location
sigma = 4 / Ts;
omega_n = sigma / zeta;
omega_d = omega_n * sqrt(1 - zeta^2);

s_desired = -8.92+ 1i *8.16 ;

disp(' ');
disp(['Desired dominant pole location: ', num2str(s_desired)]);

% 2. Calculate the Angle Deficit (using the Angle Criterion)
% Angles from all other poles and zeros, excluding the unknown compensator zero
angle_from_plant_poles = sum(angle(s_desired - plant_poles));
angle_from_plant_zeros = sum(angle(s_desired - plant_zeros));


% Sum of all known angles
known_angles_sum_rad = angle_from_plant_zeros - angle_from_plant_poles;
known_angles_sum_deg = rad2deg(known_angles_sum_rad);

% The required angle from the compensator zero to satisfy the angle criterion
% is the difference between -180 degrees and the sum of all other angles.
required_angle_from_zero_deg = -180 - known_angles_sum_deg;

% A negative angle means the zero should be to the left of the desired pole's real part.
% The angle should be positive since the desired pole is in the upper half-plane.
% We add 360 to get the correct positive angle if the result is negative.
if required_angle_from_zero_deg < 0
    required_angle_from_zero_deg = required_angle_from_zero_deg + 360;
end

required_angle_from_zero_rad = deg2rad(required_angle_from_zero_deg);

disp(['Required angle from compensator zero: ', num2str(required_angle_from_zero_deg), ' degrees']);

% 3. Calculate the Compensator Zero Location
% tan(theta) = (imaginary part) / (real part - real part of zero)
% This can be rearranged to solve for the real part of the zero.
compensator_zero_real_part = real(s_desired) - imag(s_desired) / tan(required_angle_from_zero_rad);

compensator_zero = compensator_zero_real_part;

disp(['Calculated compensator zero location: ', num2str(compensator_zero)]);