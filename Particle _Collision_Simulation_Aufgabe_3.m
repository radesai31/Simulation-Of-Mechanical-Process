%% Particle Collision Simulation via Discrete Element Method (DEM)
% This script models the interaction between two spheres using a 
% linear-hysteretic spring model and Leapfrog numerical integration.

clear; clc; close all;

% --- Configuration Parameters ---
sim_cfg.time_step = 1.0e-4;     % dt (s)
sim_cfg.duration = 2.0;         % Total time (s)
sim_cfg.restitution_coeffs = [1.0, 0.75]; % e values to evaluate

% Physical Properties
phys.mass = [0.05, 0.05];       % kg (m1, m2)
phys.diameters = [1.0, 1.0];    % m (d1, d2)
phys.k_load = 750;              % N/m (Loading stiffness)

% Initial States [x, y, z]
init.p1_pos = [0, 0, 0];
init.p2_pos = [1.1, 1.3, 0];
init.p1_vel = [0, 0, 0];
init.p2_vel = [-1.0, -1.0, 0];

% --- Execution Loop ---
for e_val = sim_cfg.restitution_coeffs
    fprintf('\n--- Analyzing Collision with Restitution e = %.2f ---\n', e_val);
    
    % Execute simulation engine
    results = run_dem_engine(sim_cfg, phys, init, e_val);
    
    [span_0](start_span)% Display numerical outputs [cite: 356-362]
    if ~isnan(results.duration)
        fprintf('  Collision Duration: %.4f s\n', results.duration);
        fprintf('  Post-Collision V1: [%.2f, %.2f, %.2f] m/s\n', results.v1_final);
        fprintf('  Post-Collision V2: [%.2f, %.2f, %.2f] m/s\n', results.v2_final);
        fprintf('  Initial Contact:   [%.2f, %.2f, %.2f] m\n', results.contact_start);
    else
        fprintf('  Warning: No contact detected.\n');
    end
end

%% --- Core Simulation Engine ---
function data = run_dem_engine(cfg, phys, init, e)
    t = 0:cfg.time_step:cfg.duration;
    n_steps = length(t);
    
    [cite_start]% Pre-allocate state matrices[span_0](end_span)
    pos1 = zeros(n_steps, 3); pos2 = zeros(n_steps, 3);
    vel1 = zeros(n_steps, 3); vel2 = zeros(n_steps, 3);
    delta = zeros(n_steps, 1);
    
    [span_1](start_span)% Set Starting Conditions[span_1](end_span)
    pos1(1,:) = init.p1_pos; pos2(1,:) = init.p2_pos;
    vel1(1,:) = init.p1_vel; vel2(1,:) = init.p2_vel;
    
    r_sum = sum(phys.diameters) / 2; % Combined radii
    
    [span_2](start_span)% Leapfrog Integration Loop[span_2](end_span)
    for i = 1:(n_steps-1)
        % 1. Calculate Forces at current step
        [F1, F2, delta(i), n_vec] = get_contact_force(pos1(i,:), pos2(i,:), vel1(i,:), vel2(i,:), phys, e, r_sum);
        
        [span_3](start_span)% 2. Predictor: Half-step velocities[span_3](end_span)
        v1_mid = vel1(i,:) + (F1 / phys.mass(1)) * (cfg.time_step / 2);
        v2_mid = vel2(i,:) + (F2 / phys.mass(2)) * (cfg.time_step / 2);
        
        [span_4](start_span)% 3. Full-step positions[span_4](end_span)
        pos1(i+1,:) = pos1(i,:) + v1_mid * cfg.time_step;
        pos2(i+1,:) = pos2(i,:) + v2_mid * cfg.time_step;
        
        [span_5](start_span)% 4. Corrector: Forces at new positions[span_5](end_span)
        [F1_next, F2_next, ~, ~] = get_contact_force(pos1(i+1,:), pos2(i+1,:), v1_mid, v2_mid, phys, e, r_sum);
        
        [span_6](start_span)% 5. Update velocities to full-step[span_6](end_span)
        vel1(i+1,:) = v1_mid + (F1_next / phys.mass(1)) * (cfg.time_step / 2);
        vel2(i+1,:) = v2_mid + (F2_next / phys.mass(2)) * (cfg.time_step / 2);
    end
    
    [span_7](start_span)% Post-processing results [cite: 411-419]
    contact_idx = find(delta > 0);
    if ~isempty(contact_idx)
        data.duration = t(contact_idx(end)) - t(contact_idx(1));
        data.v1_final = vel1(end,:);
        data.v2_final = vel2(end,:);
        
        % Geometry-based contact point
        dir_init = (pos2(contact_idx(1),:) - pos1(contact_idx(1),:));
        n_init = dir_init / norm(dir_init);
        data.contact_start = pos1(contact_idx(1),:) + (phys.diameters(1)/2) * n_init;
    else
        data.duration = NaN;
    end
    
    % Visualizations
    render_plots(t, delta, vel1, vel2, e);
end

%% --- Physics Utility Functions ---
function [f1, f2, overlap, n] = get_contact_force(p1, p2, v1, v2, phys, e, r_total)
    rel_pos = p2 - p1;
    dist = norm(rel_pos);
    [cite_start]overlap = r_total - dist;[span_7](end_span)
    
    if overlap > 0
        [span_8](start_span)n = rel_pos / dist;[span_8](end_span)
        rel_vel_n = dot(v2 - v1, n);
        
        [span_9](start_span)% Hysteretic stiffness selection [cite: 375-376]
        if rel_vel_n < 0
            k = phys.k_load; % Loading
        else
            k = (e^2) * phys.k_load; % Unloading
        end
        
        [cite_start]f_mag = max(0, k * overlap);[span_9](end_span)
        f2 = f_mag * n;
        f1 = -f2;
    else
        [span_10](start_span)[span_11](start_span)overlap = 0; n = [0,0,0]; f1 = [0,0,0]; f2 = [0,0,0]; [cite: 336-337]
    end
end

function render_plots(t, delta, v1, v2, e)
    figure('Name', sprintf('Results e=%.2f', e));
    subplot(2,1,1);
    plot(t, delta, 'k', 'LineWidth', 1.2); ylabel('Overlap (m)'); grid on;
    title(['Contact Dynamics (e = ', num2str(e), ')']);
    
    subplot(2,1,2);
    plot(t, vecnorm(v1,2,2), 'b', t, vecnorm(v2,2,2), 'r');
    xlabel('Time (s)'); ylabel('Speed (m/s)'); legend('P1','P2'); grid on;
end

