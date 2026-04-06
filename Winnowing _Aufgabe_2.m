%% Stochastic Winnowing Simulation (Monte Carlo Method)
% This script models the trajectory of grain and chaff particles through
% an air stream, accounting for randomness in particle properties.

clear; clc;

% --- Simulation Parameters ---
params.x_start = 0.50;      [span_0](start_span)% Initial x[span_0](end_span)
params.y_start = 0.50;      [span_1](start_span)% Initial y[span_1](end_span)
params.x_cutoff = 0.55;     [span_2](start_span)% Separation boundary[span_2](end_span)
params.y_floor = -0.50;     [span_3](start_span)% Ground level[span_3](end_span)
params.u_air = 0.2;         [span_4](start_span)% Reference air velocity[span_4](end_span)
params.h_nozzle = 0.10;     [span_5](start_span)% Nozzle height[span_5](end_span)

% Physical Constants
env.rho_air = 1.2;          [span_6](start_span)% Air density[span_6](end_span)
env.mu_air = 1.8e-5;        [span_7](start_span)% Air viscosity[span_7](end_span)
env.g = 9.81;               [span_8](start_span)% Gravity[span_8](end_span)

% Monte Carlo Settings
mc.n_samples = 1000;        [span_9](start_span)% Number of particles to simulate[span_9](end_span)
mc.n_internal = 10;         [span_10](start_span)% Internal MC steps for time averaging[span_10](end_span)
mc.dt = 0.005;              [span_11](start_span)[span_12](start_span)% Time step[span_11](end_span)[span_12](end_span)

% --- Particle Distribution Functions ---
% Chaff: Uniform diameter, Normal density
get_chaff_dia = @(n) 0.002 + (0.005 - 0.002) * rand(1, n); [span_13](start_span)%
get_chaff_rho = @(n) 50 + 20 * randn(1, n);                %[span_13](end_span)

% Grain: Normal diameter, Constant density
get_grain_dia = @(n) 0.0025 + 0.001 * randn(1, n);         [span_14](start_span)%
rho_grain_fixed = 750;                                     %[span_14](end_span)

% --- Physics Functions ---
[span_15](start_span)% Local fluid velocity profile[span_15](end_span)
calc_v_fluid = @(px, py) 6.2 * params.u_air * sqrt(params.h_nozzle / px) * ...
                         exp(-50 * (py^2 / px^2));

[span_16](start_span)% Particle Terminal Velocity[span_16](end_span)
calc_v_term = @(p_rho, p_dia) sqrt((4/3) * (p_rho - env.rho_air) * env.g * ...
                                   p_dia / (3 * env.rho_air));

% --- Execution Logic ---
results.bin_left = 0;
results.bin_right = 0;

for type = 1:2 % 1 = Grain, 2 = Chaff
    % Generate distributions for the batch
    if type == 1
        d_list = get_grain_dia(mc.n_samples);
        r_list = repmat(rho_grain_fixed, 1, mc.n_samples);
    else
        d_list = get_chaff_dia(mc.n_samples);
        r_list = get_chaff_rho(mc.n_samples);
    end
    
    [span_17](start_span)% Random launch angles (-95 to -85 degrees)[span_17](end_span)
    launch_angles = -95 + (10) * rand(1, mc.n_samples);
    
    for m = 1:mc.n_samples
        % Particle-specific properties
        d = d_list(m);
        rho = r_list(m);
        vol = (pi * d^3) / 6; [span_18](start_span)%
        
        % Initial State based on terminal velocity[span_18](end_span)
        vt = calc_v_term(rho, d);
        pos = [params.x_start, params.y_start];
        vel = [abs(vt * cosd(180 - launch_angles(m))), ...
               abs(vt * cosd(launch_angles(m) + 90))];
        
        % Path Integration
        [span_19](start_span)[span_20](start_span)while pos(2) > params.y_floor %[span_19](end_span)[span_20](end_span)
            % Local environment
            vf = calc_v_fluid(pos(1), pos(2));
            v_rel = [vf - vel(1), -vel(2)];
            v_mag = norm(v_rel); [span_21](start_span)%
            
            % Drag Coefficient[span_21](end_span)
            Re = (env.rho_air * v_mag * d) / env.mu_air;
            Cd = (Re < 800) * (24/Re * (1 + 0.15*Re^0.687)) + (Re >= 800) * 0.44;
            
            [span_22](start_span)[span_23](start_span)% Forces[span_22](end_span)[span_23](end_span)
            f_drag = 0.5 * pi * d^2 * env.rho_air * Cd * v_mag * v_rel;
            f_grav = rho * vol * env.g;
            f_buoy = -env.rho_air * vol * env.g;
            
            [span_24](start_span)[span_25](start_span)% Monte Carlo Time Averaging[span_24](end_span)[span_25](end_span)
            mc_sum = sum(0 + rand(1, mc.n_internal) * mc.dt); 
            mc_avg_t = mc_sum / mc.n_internal;
            
            [span_26](start_span)% Euler Update [cite: 65-68, 108-111]
            pos(1) = pos(1) + mc.dt * vel(1) * mc_avg_t;
            pos(2) = pos(2) - mc.dt * vel(2) * mc_avg_t;
            
            accel_x = f_drag(1) / (rho * vol);
            accel_y = (f_grav + f_buoy + f_drag(2)) / (rho * vol);
            
            vel(1) = vel(1) + mc.dt * accel_x * mc_avg_t;
            vel(2) = vel(2) + mc.dt * accel_y * mc_avg_t;
        end
        
        [cite_start]% Sort into bins[span_26](end_span)
        if pos(1) < params.x_cutoff
            results.bin_left = results.bin_left + 1;
        else
            results.bin_right = results.bin_right + 1;
        end
    end
end

% --- Final Statistics ---
fprintf('Particle Distribution Analysis:\n');
fprintf('Bin 1 (Left): %.2f%%\n', (results.bin_left / (2 * mc.n_samples)) * 100);
fprintf('Bin 2 (Right): %.2f%%\n', (results.bin_right / (2 * mc.n_samples)) * 100);

