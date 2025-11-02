%% SIMULINK MODEL SETUP FOR SELF-BALANCING ROBOT
% This script creates a Simulink model and provides guidance on setup
% Run this after robot_model.m to create the complete Simulink system

%% ========================================================================
% SETUP: Create Simulink Model Components
% =========================================================================

% Create a new Simulink model
model_name = 'robot_balance_system';
if bdIsLoaded(model_name)
    close_system(model_name, 0);
end

new_system(model_name);
load_system(model_name);

% Define block positions [x, y, x+width, y+height]
pos_in = [50, 150, 100, 170];
pos_ctrl = [150, 100, 250, 200];
pos_plant = [300, 100, 450, 200];
pos_out = [550, 150, 600, 170];
pos_scope = [550, 50, 650, 150];
pos_mux = [450, 100, 500, 200];

%% Add Simulink Blocks for Linear System

% Input: Step torque command
add_block('simulink/Sources/Step', [model_name '/Torque Command']);
set_param([model_name '/Torque Command'], 'Time', '1', 'After', '5', 'Position', pos_in);

% Gain block for motor torque scaling (if needed)
add_block('simulink/Math Operations/Gain', [model_name '/Torque Gain']);
set_param([model_name '/Torque Gain'], 'Gain', '10', 'Position', pos_ctrl);

% State-Space block for plant dynamics
add_block('simulink/Continuous/State-Space', [model_name '/Plant (Linear)']);
% Note: After running robot_model.m, A, B, C, D are available in workspace
set_param([model_name '/Plant (Linear)'], 'Position', pos_plant);
set_param([model_name '/Plant (Linear)'], 'A', 'A', 'B', 'B', 'C', 'C', 'D', 'D');
set_param([model_name '/Plant (Linear)'], 'InitialConditionSource', 'external');
set_param([model_name '/Plant (Linear)'], 'X0', '[0; 0; 0.1; 0]');

% Mux to combine state outputs
add_block('simulink/Signal Routing/Mux', [model_name '/Mux']);
set_param([model_name '/Mux'], 'Inputs', '4', 'Position', pos_mux);

% Scope for monitoring outputs
add_block('simulink/Sinks/Scope', [model_name '/State Monitor']);
set_param([model_name '/State Monitor'], 'Position', pos_scope);
set_param([model_name '/State Monitor'], 'NumInputPorts', '4');

% Add output display
add_block('simulink/Sinks/To Workspace', [model_name '/State Output']);
set_param([model_name '/State Output'], 'VariableName', 'z_out');
set_param([model_name '/State Output'], 'Position', [550, 200, 650, 220]);

%% Connect Blocks with Lines

% Torque Command -> Gain
add_line(model_name, 'Torque Command/1', 'Torque Gain/1');

% Gain -> Plant
add_line(model_name, 'Torque Gain/1', 'Plant (Linear)/1');

% Plant outputs -> Mux
add_line(model_name, 'Plant (Linear)/1', 'Mux/1');

% Mux -> Scope and To Workspace
add_line(model_name, 'Mux/1', 'State Monitor/1');
add_line(model_name, 'Mux/1', 'State Output/1');

% Auto-route remaining connections
set_param(model_name, 'EnableLBReport', 'on');
set_param(model_name, 'ShowPortDataTypes', 'on');

fprintf('Simulink model created: %s\n', model_name);
fprintf('Block connections established\n\n');

%% ========================================================================
% PART 2: CONFIGURE SOLVER SETTINGS
% =========================================================================

set_param(model_name, 'Solver', 'ode45');
set_param(model_name, 'StartTime', '0');
set_param(model_name, 'StopTime', '10');
set_param(model_name, 'SolverType', 'Variable-step');
set_param(model_name, 'MaxStep', '0.01');

fprintf('Solver configured:\n');
fprintf('  Method: ode45 (Runge-Kutta 4-5)\n');
fprintf('  Stop Time: 10 seconds\n');
fprintf('  Max Step: 0.01 seconds\n\n');

%% ========================================================================
% PART 3: ALTERNATIVE - NONLINEAR S-FUNCTION MODEL
% =========================================================================

% Create S-function for nonlinear dynamics
s_func_code = [
    '%% S-Function for Robot Nonlinear Dynamics\n'
    'function [sys,x0,str,ts,xts] = robot_sfun(t,x,u,flag,params)\n'
    '\n'
    'switch flag\n'
    '  case 0    % Initialization\n'
    '    [sys,x0,str,ts,xts] = mdlInitializeSizes(params);\n'
    '\n'
    '  case 1    % Derivatives\n'
    '    sys = mdlDerivatives(t,x,u,params);\n'
    '\n'
    '  case 3    % Outputs\n'
    '    sys = mdlOutputs(t,x,u);\n'
    '\n'
    '  case {2,4,9}\n'
    '    sys = [];\n'
    '\n'
    '  otherwise\n'
    '    DAStudio.error(''Simulink:InvalidFlag'',num2str(flag));\n'
    'end\n'
    '\n'
    'function [sys,x0,str,ts,xts] = mdlInitializeSizes(params)\n'
    '  sys = simsizes;\n'
    '  sys.NumContStates  = 4;      % 4 continuous states\n'
    '  sys.NumDiscStates  = 0;\n'
    '  sys.NumOutputs     = 4;      % Output all states\n'
    '  sys.NumInputs      = 1;      % Motor torque input\n'
    '  sys.DirFeedthrough = 0;\n'
    '  sys.NumSampleTimes = 1;\n'
    '  x0 = [0; 0; 0.1; 0];         % Initial state\n'
    '  str = [];\n'
    '  ts = [0 0];\n'
    '  xts = [];\n'
    '\n'
    'function sys = mdlDerivatives(t,x,u,params)\n'
    '  % Extract state\n'
    '  x_pos = x(1);\n'
    '  x_dot = x(2);\n'
    '  theta = x(3);\n'
    '  theta_dot = x(4);\n'
    '  \n'
    '  % Extract parameters\n'
    '  alpha = params.alpha;\n'
    '  beta = params.beta;\n'
    '  gamma = params.gamma;\n'
    '  g = params.g;\n'
    '  R = params.R;\n'
    '  b = params.b;\n'
    '  \n'
    '  % Nonlinear dynamics\n'
    '  Delta = alpha*gamma - beta^2*cos(theta)^2;\n'
    '  \n'
    '  x_ddot = (1/Delta) * (gamma*(u/R + beta*theta_dot^2*sin(theta)) - ...\n'
    '           beta*cos(theta)*(-u + g*beta*sin(theta) - b*theta_dot));\n'
    '  \n'
    '  theta_ddot = (1/Delta) * (alpha*(-u + g*beta*sin(theta) - b*theta_dot) - ...\n'
    '                beta*cos(theta)*(u/R + beta*theta_dot^2*sin(theta)));\n'
    '  \n'
    '  sys = [x_dot; x_ddot; theta_dot; theta_ddot];\n'
    '\n'
    'function sys = mdlOutputs(t,x,u)\n'
    '  sys = x;  % Output all states\n'
];

% Save S-function to file
fid = fopen('robot_sfun.m', 'w');
fprintf(fid, '%s', s_func_code);
fclose(fid);
fprintf('S-function file created: robot_sfun.m\n\n');

%% ========================================================================
% PART 4: SETUP INSTRUCTIONS FOR MANUAL SIMULINK BUILDING
% =========================================================================

fprintf('=====================================================\n');
fprintf('MANUAL SIMULINK MODEL SETUP INSTRUCTIONS\n');
fprintf('=====================================================\n\n');

fprintf('To create a complete Simulink model manually:\n\n');

fprintf('STEP 1: INPUT BLOCKS\n');
fprintf('-  Step Input: Set Time=1, After=5, Initial=0, Final=5\n');
fprintf('   (This creates a step torque command starting at t=1s)\n\n');

fprintf('STEP 2: CONTROLLER BLOCK (Optional for testing)\n');
fprintf('   A. PID Controller Block:\n');
fprintf('      - Set proportional, integral, derivative gains\n');
fprintf('      - Suggested: Kp=50, Ki=10, Kd=5 (tune from here)\n\n');
fprintf('   B. Or State Feedback Controller:\n');
fprintf('      - Multiply states by feedback gain matrix K\n');
fprintf('      - K can be computed using LQR design\n\n');

fprintf('STEP 3: PLANT MODEL\n');
fprintf('   A. Linear Model (for control design):\n');
fprintf('      - State-Space block with A, B, C, D matrices\n');
fprintf('      - Initial condition: [0; 0; 0.1; 0]\n\n');
fprintf('   B. Nonlinear Model (for validation):\n');
fprintf('      - S-Function block pointing to robot_sfun.m\n');
fprintf('      - Same initial conditions as linear model\n\n');

fprintf('STEP 4: OUTPUT AND MONITORING\n');
fprintf('   - Scope blocks to monitor: angle, angular velocity\n');
fprintf('   - To Workspace blocks to save data for analysis\n');
fprintf('   - Display blocks to show real-time values\n\n');

fprintf('STEP 5: SOLVER SETTINGS\n');
fprintf('   - Go to Simulation > Model Configuration Parameters\n');
fprintf('   - Solver: ode45 (or ode23t for stiff systems)\n');
fprintf('   - Max Step Size: 0.001 to 0.01 seconds\n');
fprintf('   - Stop Time: 10-15 seconds\n\n');

fprintf('STEP 6: RUN SIMULATION\n');
fprintf('   - Execute: sim(model_name) or press Run button\n');
fprintf('   - Plot and analyze results\n\n');

%% ========================================================================
% PART 5: LQR CONTROLLER DESIGN (for simulation)
% =========================================================================

fprintf('=====================================================\n');
fprintf('LQR CONTROLLER DESIGN\n');
fprintf('=====================================================\n\n');

% Define LQR cost matrices
Q_lqr = diag([1, 1, 1000, 1]);  % Higher weight on angle error
R_lqr = 1;                       % Lower control effort cost

fprintf('LQR Design:\n');
fprintf('  Q matrix (state cost): diag([1, 1, 1000, 1])\n');
fprintf('  R matrix (input cost): 1\n');
fprintf('  Higher Q on angle to prioritize balance\n\n');

% Note: K can be computed as:
% [K, S, e] = lqr(A, B, Q_lqr, R_lqr);
% fprintf('Feedback gain K (computed offline):\n');
% disp(K);

fprintf('To implement LQR controller in Simulink:\n');
fprintf('1. Compute K = lqr(A, B, Q, R) in MATLAB workspace\n');
fprintf('2. Add Gain block with value -K\n');
fprintf('3. Multiply state vector by -K to get control torque\n');
fprintf('4. Command: tau = -K * z, where z is state vector\n\n');

%% ========================================================================
% PART 6: EXAMPLE TESTS TO RUN
% =========================================================================

fprintf('=====================================================\n');
fprintf('RECOMMENDED SIMULINK TESTS\n');
fprintf('=====================================================\n\n');

fprintf('TEST 1: Open-Loop Response\n');
fprintf('  - No controller, just step torque input\n');
fprintf('  - Observe system instability (angle diverges)\n\n');

fprintf('TEST 2: Impulse Disturbance Rejection\n');
fprintf('  - Apply brief torque pulse (0-0.1s)\n');
fprintf('  - With controller, system should reject disturbance\n\n');

fprintf('TEST 3: Step Angle Reference Tracking\n');
fprintf('  - Command desired angle via controller\n');
fprintf('  - Observe settling time and overshoot\n\n');

fprintf('TEST 4: Sinusoidal Reference Tracking\n');
fprintf('  - Frequency sweep from 0.1 to 10 Hz\n');
fprintf('  - Compare amplitude and phase response\n\n');

fprintf('TEST 5: Nonlinear vs Linear Model Comparison\n');
fprintf('  - Run same input on both models\n');
fprintf('  - Verify linearization validity for small angles\n\n');

fprintf('TEST 6: Robustness to Parameter Variations\n');
fprintf('  - Vary mass m_b by ±20%%\n');
fprintf('  - Vary height h by ±10%%\n');
fprintf('  - Test controller stability margins\n\n');

%% ========================================================================
% PART 7: VISUALIZATION FUNCTION
% =========================================================================

fprintf('Creating visualization helper function...\n');
fid = fopen('plot_robot_states.m', 'w');
fprintf(fid, [
    'function plot_robot_states(t, z, u)\n'
    '%% PLOT_ROBOT_STATES Visualize simulation results\n'
    '%  t - time vector\n'
    '%  z - state matrix [x, x_dot, theta, theta_dot]\n'
    '%  u - input vector (torque)\n'
    '\n'
    'figure(''Name'', ''Robot States'', ''NumberTitle'', ''off'');\n'
    '\n'
    'subplot(4,1,1);\n'
    'plot(t, u, ''LineWidth'', 2); grid on;\n'
    'ylabel(''Motor Torque (N⋅m)'');\n'
    'title(''Self-Balancing Robot Simulation Results'');\n'
    '\n'
    'subplot(4,1,2);\n'
    'plot(t, z(:,1), ''LineWidth'', 2); grid on;\n'
    'ylabel(''Position x (m)'');\n'
    '\n'
    'subplot(4,1,3);\n'
    'plot(t, z(:,3)*180/pi, ''LineWidth'', 2); grid on;\n'
    'ylabel(''Angle θ (deg)'');\n'
    'axhline(0, ''Color'', ''r'', ''LineStyle'', ''--'', ''Alpha'', 0.5);\n'
    '\n'
    'subplot(4,1,4);\n'
    'plot(t, z(:,2), ''LineWidth'', 2); hold on;\n'
    'plot(t, z(:,4)*180/pi, ''LineWidth'', 2); grid on;\n'
    'ylabel(''Velocities'');\n'
    'xlabel(''Time (s)'');\n'
    'legend(''ẋ (m/s)'', ''ω (deg/s)'');\n'
    '\n'
    'end\n'
]);
fclose(fid);
fprintf('Visualization function created: plot_robot_states.m\n\n');

fprintf('=====================================================\n');
fprintf('SETUP COMPLETE!\n');
fprintf('=====================================================\n');
fprintf('Files created:\n');
fprintf('  1. Simulink model: %s.slx (auto-generated)\n', model_name);
fprintf('  2. S-function: robot_sfun.m\n');
fprintf('  3. Plotting tool: plot_robot_states.m\n\n');

fprintf('Next steps:\n');
fprintf('1. Run: robot_model.m (if not already done)\n');
fprintf('2. Open Simulink model: %s.slx\n', model_name);
fprintf('3. Configure simulation parameters\n');
fprintf('4. Run simulations with various test inputs\n');
fprintf('5. Analyze results using plot_robot_states.m\n\n');
