% --- Load the Bode plot figure data ---
fig_handle = openfig('bode_q1.fig');
all_axes = findall(fig_handle, 'type', 'axes'); 

% --- Locate the specific axis for "Magnitude (dB)" ---
magnitude_axis = []; 
for i = 1:length(all_axes)
    current_ylabel = get(get(all_axes(i), 'YLabel'), 'String');
    if contains(lower(current_ylabel), 'mag')
        magnitude_axis = all_axes(i); 
        break;
    end
end

% --- Extract frequency and magnitude data from the plot line ---
line_objects = findall(magnitude_axis, 'type', 'line');
data_line = line_objects(1); 
fig_frequencies  = get(data_line, 'XData'); 
fig_magnitudes_db = get(data_line, 'YData'); 
close(fig_handle); 

% --- Define the data range to be used for fitting ---

frequencies_for_fit = fig_frequencies; 
magnitudes_for_fit  = fig_magnitudes_db;  

% --- Define search space for model parameters ---
best_mse = inf; 
best_parameters = [NaN NaN NaN]; 

% Define the parameter ranges for the brute-force sweep
wn_range   = 4.3:0.05:4.7;    
K_range    = 0.5:0.05:0.8;    
zeta_range = 0.1:0.001:0.12;    

% --- Initialize plot for visualization ---
figure; 
hold on;
grid on;

% --- Brute-force parameter sweep ---
fprintf('Starting parameter sweep...\n');
plot_color = [0.85 0.85 0.85]; 
for wn = wn_range
    for zeta = zeta_range
        for K = K_range
            % --- Define the system transfer function ---
            % Model: G(s) = K * wn^2 / (s * (s^2 + 2*zeta*wn*s + wn^2))
            num = K * wn^2;
            den = [1, 2*zeta*wn, wn^2, 0];
            current_model = tf(num, den);

            % --- Calculate model's magnitude response ---
            % We only need one calculation for both MSE and plotting
            [model_mag_linear, ~] = bode(current_model, frequencies_for_fit);
            model_mag_db = 20*log10(squeeze(model_mag_linear));
            
            % --- Compute Mean Squared Error (MSE) ---
            squared_errors = (model_mag_db - magnitudes_for_fit).^2;
            current_mse = mean(squared_errors(~isnan(squared_errors))); % Ignore NaNs

            % --- Plot the current sweep's response ---
            % (This plots every iteration for a "ghosted" effect)
            plot(fig_frequencies, model_mag_db, 'Color', plot_color);

            % --- Check if this model is the new best fit ---
            if current_mse < best_mse
                best_mse = current_mse;
                best_parameters = [K, zeta, wn];
            end
        end
    end
end
fprintf('Parameter sweep finished.\n');

% --- Extract best parameters from the search ---
K_best = best_parameters(1);
zeta_best = best_parameters(2);
wn_best = best_parameters(3);

% --- Display the best-fit results ---
fprintf('Best fit parameters:\n');
fprintf('  K     = %.4f\n', K_best);
fprintf('  zeta  = %.4f\n', zeta_best);
fprintf('  wn    = %.4f rad/s\n', wn_best);
fprintf('  MSE   = %.6f\n', best_mse);

% --- Generate the final best-fit model response ---
best_fit_model = tf(K_best * wn_best^2, [1, 2*zeta_best*wn_best, wn_best^2, 0]);
[best_model_mag_linear, ~] = bode(best_fit_model, fig_frequencies);
best_model_mag_db = 20*log10(squeeze(best_model_mag_linear));

% --- Finalize plot: Overlay data and best-fit model ---
sweep_plot_handles = findall(gca, 'Color', plot_color); 
measured_data_plot = plot(frequencies_for_fit, magnitudes_for_fit, 'b', 'LineWidth', 2.5); 
best_fit_plot = plot(fig_frequencies, best_model_mag_db, 'r--', 'LineWidth', 2); 

% --- Apply final formatting to the plot ---
set(gca, 'XScale', 'log'); 
xlim([1e-1 1e2]);          
ylim([-60, 40]);           
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Bode Magnitude: Model Fit vs Measured Data');
legend([sweep_plot_handles(1), measured_data_plot, best_fit_plot], ...
       'Model sweeps', 'Measured data', 'Best-fit model', 'Location', 'southwest');
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
hold off;