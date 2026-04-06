%% Particle Trajectory Simulation (Winnowing Process)
% This script simulates the movement of grain and chaff in an air stream
% using the Runge-Kutta (RK4) numerical integration method.

clear; clc;

% --- Configuration & Constants ---
start_pos = [0.50, 0.50];     % Initial [x, y] coordinates
bin_threshold = -0.50;        % Y-coordinate where bins are located
u_ref = 0.1;                  % Reference fluid velocity
nozzle_h = 0.10;              % Nozzle height parameter

% Environmental Properties
air_density = 1.2;            % rho_f (kg/m^3)
air_viscosity = 1.8e-5;       % mu_f (Pa*s)
gravity = 9.81;

% Particle Definitions [Density, Diameter]
grain_props = [750, 2.5e-3]; 
chaff_props = [50, 3.25e-3];

% Simulation Parameters
num_steps = 200;
time_increment = 1 / num_steps;

% --- Data Storage Setup ---
% Results stored as: trajectories{1} = Grain, trajectories{2} = Chaff
trajectories = cell(1, 2); 

% --- Physics Model Functions ---
% Fluid velocity profile (x-direction only)
get_v_fluid = @(px, py) 6.2 * u_ref * sqrt(nozzle_h / px) * exp(-50 * (py^2 / px^2));

% --- Main Simulation Loop ---
for p_type = 1:2
    % Select properties based on iteration
    if p_type == 1
        p_rho = grain_props(1); p_dia = grain_props(2);
    else
        p_rho = chaff_props(1); p_dia = chaff_props(2);
    end
    
    % Derived Geometric properties
    p_vol = (pi * p_dia^3) / 6;
    p_area = pi * (p_dia/2)^2;
    
    % Initial States: [x, y, vx, vy]
    state = [start_pos(1), start_pos(2), 0, 0];
    history = state;
    t = 0;
    
    % Time Integration (Falling until bin threshold reached)
    while state(2) > bin_threshold
        curr_x = state(1);
        curr_y = state(2);
        curr_vx = state(3);
        curr_vy = state(4);
        
        % Calculate Forces at current step
        v_f = get_v_fluid(curr_x, curr_y);
        rel_v_x = v_f - curr_vx;
        rel_v_y = 0 - curr_vy;
        v_mag = sqrt(rel_v_x^2 + rel_v_y^2);
        
        % Reynolds number and Drag Coefficient
        Re = (air_density * v_mag * p_dia) / air_viscosity;
        if Re < 800
            coeff_d = (24/Re) * (1 + 0.15 * Re^0.687);
        else
            coeff_d = 0.44;
        end
        
        % Force Calculations
        f_drag_x = 0.5 * air_density * p_area * coeff_d * v_mag * rel_v_x;
        f_drag_y = 0.5 * air_density * p_area * coeff_d * v_mag * rel_v_y;
        f_net_y = (p_rho * p_vol * gravity) - (air_density * p_vol * gravity) + f_drag_y;
        
        % Define Derivatives for RK4: d[x, y, vx, vy]/dt
        % Note: y-position decreases, so we treat gravity-driven motion accordingly
        derivs = @(s) [s(3), -s(4), f_drag_x/(p_rho*p_vol), f_net_y/(p_rho*p_vol)];
        
        % RK4 Integration Steps
        k1 = derivs(state);
        k2 = derivs(state + 0.5 * time_increment * k1);
        k3 = derivs(state + 0.5 * time_increment * k2);
        k4 = derivs(state + time_increment * k3);
        
        % Update State
        state = state + (time_increment / 6) * (k1 + 2*k2 + 2*k3 + k4);
        history = [history; state]; 
    end
    trajectories{p_type} = history;
end

% --- Visualization ---
figure('Color', 'w');
plot(trajectories{1}(:,1), trajectories{1}(:,2), 'LineWidth', 2, 'Color', [0 0.4 0.8]);
hold on;
plot(trajectories{2}(:,1), trajectories{2}(:,2), 'LineWidth', 2, 'Color', [0.8 0.2 0]);
grid on;

% Styling and Labels
xlabel('Horizontal Distance (x)');
ylabel('Vertical Height (y)');
title('Comparative Particle Trajectories in Winnowing Flow');
legend('Grain', 'Chaff', 'Location', 'best');
set(gca, 'FontSize', 10);

