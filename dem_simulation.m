function [overlap, normal_vector] = calculate_overlap_and_normal(pos1, diameter1, pos2, diameter2)
%CALCULATE_OVERLAP_AND_NORMAL Computes the overlap and normal vector between two particles.
%
%   [overlap, normal_vector] = calculate_overlap_and_normal(pos1, diameter1, pos2, diameter2)
%
%   Inputs:
%     pos1       : 3D position vector of the center of particle 1 (e.g., [x1, y1, z1])
%     diameter1  : Diameter of particle 1 (scalar)
%     pos2       : 3D position vector of the center of particle 2 (e.g., [x2, y2, z2])
%     diameter2  : Diameter of particle 2 (scalar)
%
%   Outputs:
%     overlap        : Scalar value representing the penetration depth between particles.
%                      Returns 0 if particles are not overlapping.
%     normal_vector  : 3D unit vector pointing from particle 1 to particle 2 along
%                      the line of centers. Returns [0,0,0] if no overlap.

    % Calculate radii
    radius1 = diameter1 / 2;
    radius2 = diameter2 / 2;

    % Desired distance between centers for contact (sum of radii)
    desired_distance = radius1 + radius2;

    % Vector connecting the centers of the two particles
    r_vec = pos2 - pos1;

    % Current distance between centers
    current_distance = norm(r_vec); % norm() calculates the Euclidean distance (magnitude)

    % Calculate overlap
    overlap = desired_distance - current_distance;

    % Determine the normal vector
    if overlap > 0
        % Particles are overlapping, calculate unit normal vector
        normal_vector = r_vec / current_distance;
    else
        % No overlap, set overlap to 0 and normal vector to zero vector
        overlap = 0;
        normal_vector = [0, 0, 0];
    end

end


    clear; clc; close all; % Clear workspace, command window, close figures

    fprintf('Starting combined DEM simulation and analysis...\n');

    % --- Define Common Simulation Parameters ---
    dt = 1.0e-4; % s, Time step (small enough for stable integration)
    t_end = 2.0; % s, End time of simulation (adjust as needed to capture full collision)
    
    % --- Define Common Particle Properties ---
    m1 = 0.05; % kg, Mass of particle 1
    m2 = 0.05; % kg, Mass of particle 2
    d1 = 1.0;  % m, Diameter of particle 1 (Radius R1 = 0.5 m)
    d2 = 1.0;  % m, Diameter of particle 2 (Radius R2 = 0.5 m)
    
    % Material property for DEM (from problem)
    k_loading = 750; % N/m, Normal loading spring constant (k_loading_n)

    % --- Define Common Initial Conditions ---
    x1_initial = [0, 0, 0]; % m
    x2_initial = [1.1, 1.3, 0]; % m
    v1_initial = [0, 0, 0]; % m/s
    v2_initial = [-1.0, -1.0, 0]; % m/s

    % --- Define the 'e' values to test for Problem 5 ---
    e_values_to_test = [1.0, 0.75]; % Run simulation for e=1.0 (P5b) and e=0.75 (P5c/d)

    % =====================================================================
    % --- PROBLEM 5 CODE: Orchestration and Display of Results ---
    % This section controls running the simulation for different 'e' values
    % and prints the specific results requested in Problem 5.
    % =====================================================================

    for k = 1:length(e_values_to_test)
        current_e_for_dem = e_values_to_test(k);

        fprintf('\n=================================================================\n');
        fprintf('                     RUNNING SIMULATION for e = %.2f                \n', current_e_for_dem);
        fprintf('=================================================================\n');

        % Call the core simulation logic (Problem 4, defined below as a nested function)
        [v1_post, v2_post, collision_duration, initial_contact_pt, final_contact_pt] = ...
            runCoreDemSimulation(dt, t_end, m1, m2, d1, d2, k_loading, current_e_for_dem, ...
                                 x1_initial, x2_initial, v1_initial, v2_initial);

        % --- Display Results for Problem 5 based on current 'e' ---
        fprintf('\n--- Numerical Results for Problem 5 (for e = %.2f) ---\n', current_e_for_dem);
        
        if current_e_for_dem == 1.0
            fprintf('  (Results for Problem 5b - e=1.0)\n');
        elseif current_e_for_dem == 0.75
            fprintf('  (Results for Problem 5c, 5d - e=0.75)\n');
        end

        fprintf('  Post-collision velocity P1: [%f, %f, %f] m/s\n', v1_post(1), v1_post(2), v1_post(3));
        fprintf('  Post-collision velocity P2: [%f, %f, %f] m/s\n', v2_post(1), v2_post(2), v2_post(3));
        
        if ~isnan(collision_duration)
            fprintf('  Collision Duration: %f seconds\n', collision_duration); % Problem 5d
        else
            fprintf('  Collision Duration: Not applicable (no overlap detected)\n');
        end

        fprintf('  Initial contact point: [%f, %f, %f] m (Problem 5a conceptual)\n', ...
            initial_contact_pt(1), initial_contact_pt(2), initial_contact_pt(3));
        fprintf('  Contact point at separation: [%f, %f, %f] m (Problem 5e conceptual)\n', ...
            final_contact_pt(1), final_contact_pt(2), final_contact_pt(3));

    end

    fprintf('\nAll simulations and analyses completed.\n');

    % =====================================================================
    % --- PROBLEM 4 CODE: Core DEM Simulation Logic ---
    % This nested function contains the main simulation loop, implementing
    % the Leapfrog algorithm and the hysteretic contact force model.
    % =====================================================================

    function [v1_post_collision, v2_post_collision, duration_of_collision, contact_point_initial, contact_point_final] = ...
        runCoreDemSimulation(dt, t_end, m1, m2, d1, d2, k_loading, e_for_dem, ...
                             x1_initial, x2_initial, v1_initial, v2_initial)

        % Re-calculate values that depend on common parameters but are used within the simulation
        r1 = d1/2; % m, Radius of particle 1
        r2 = d2/2; % m, Radius of particle 2
        t_vec = 0:dt:t_end; % Time vector for this specific run
        num_steps = length(t_vec);

        % Initialize Arrays for Storing Results
        x1_history = zeros(num_steps, 3);
        x2_history = zeros(num_steps, 3);
        v1_history = zeros(num_steps, 3);
        v2_history = zeros(num_steps, 3);
        overlap_history = zeros(num_steps, 1);

        % Set initial conditions
        x1_history(1,:) = x1_initial;
        x2_history(1,:) = x2_initial;
        v1_history(1,:) = v1_initial;
        v2_history(1,:) = v2_initial;

        current_x1 = x1_initial;
        current_x2 = x2_initial;
        current_v1 = v1_initial;
        current_v2 = v2_initial;

        % Main Simulation Loop (Leapfrog Integration)
        for i = 1:num_steps-1

            % Calculate Forces at current time t
            % Calls the nested function for Problem 3
            [overlap, n_hat] = calculate_overlap_and_normal(current_x1, d1, current_x2, d2);
            overlap_history(i) = overlap; % Store overlap history

            F1_total = [0,0,0];
            F2_total = [0,0,0];

            if overlap > 0 % Check if particles are in contact (overlapping)
                v_rel = current_v2 - current_v1;
                v_rel_n = dot(v_rel, n_hat);

                % Determine current stiffness based on loading/unloading phase
                if v_rel_n < 0 % Loading phase: particles approaching (overlap increasing)
                    current_k_contact = k_loading;
                else % Unloading phase: particles separating (overlap decreasing) or at peak compression
                    current_k_contact = (e_for_dem^2) * k_loading;
                end

                F_n_magnitude = current_k_contact * overlap;
                if F_n_magnitude < 0 % Ensure force is always repulsive
                    F_n_magnitude = 0;
                end
                F_normal_vec = F_n_magnitude * n_hat;

                F1_total = -F_normal_vec;
                F2_total = F_normal_vec;
            end

            % Calculate Accelerations at current time t
            a1 = F1_total / m1;
            a2 = F2_total / m2;

            % Leapfrog Predictor: Update velocities to half-step (v(t+dt/2))
            v1_half_step = current_v1 + a1 * (dt / 2);
            v2_half_step = current_v2 + a2 * (dt / 2);

            % Update positions to full step (x(t+dt)) using half-step velocities
            next_x1 = current_x1 + v1_half_step * dt;
            next_x2 = current_x2 + v2_half_step * dt;
            
            % Store predicted positions
            x1_history(i+1,:) = next_x1;
            x2_history(i+1,:) = next_x2;

            % Re-calculate Forces at next predicted position (t+dt) for corrector
            % Calls the nested function for Problem 3
            [overlap_next, n_hat_next] = calculate_overlap_and_normal(next_x1, d1, next_x2, d2);

            F1_total_next = [0,0,0];
            F2_total_next = [0,0,0];

            if overlap_next > 0
                v_rel_next_pred = v2_half_step - v1_half_step;
                v_rel_n_next = dot(v_rel_next_pred, n_hat_next);

                if v_rel_n_next < 0
                    current_k_contact_next = k_loading;
                else
                    current_k_contact_next = (e_for_dem^2) * k_loading;
                end

                F_n_magnitude_next = current_k_contact_next * overlap_next;
                if F_n_magnitude_next < 0
                    F_n_magnitude_next = 0;
                end
                F_normal_vec_next = F_n_magnitude_next * n_hat_next;

                F1_total_next = -F_normal_vec_next;
                F2_total_next = F_normal_vec_next;
            end

            % Calculate accelerations at next predicted time (t+dt)
            a1_next = F1_total_next / m1;
            a2_next = F2_total_next / m2;

            % Leapfrog Corrector: Update velocities to full step (v(t+dt))
            next_v1 = v1_half_step + a1_next * (dt / 2);
            next_v2 = v2_half_step + a2_next * (dt / 2);

            % Store updated velocities
            v1_history(i+1,:) = next_v1;
            v2_history(i+1,:) = next_v2;

            % Update current values for next iteration
            current_x1 = next_x1;
            current_x2 = next_x2;
            current_v1 = next_v1;
            current_v2 = next_v2;

            % Display progress (optional)
            if mod(i, round(num_steps/10)) == 0
                fprintf('  Progress for e=%.2f: %.1f%%\n', e_for_dem, (i/num_steps)*100);
            end
        end
        fprintf('  Simulation for e=%.2f finished.\n', e_for_dem);

        % --- Plotting for the current 'e' value ---
        t_vec_full = 0:dt:t_end; % Re-create time vector for plotting

        % Plot Particle Trajectories (X-Y plane)
        figure('Name', sprintf('Particle Trajectories (e = %.2f)', e_for_dem));
        plot(x1_history(:,1), x1_history(:,2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Particle 1 Path'); hold on;
        plot(x2_history(:,1), x2_history(:,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Particle 2 Path');
        % Mark initial positions
        plot(x1_initial(1), x1_initial(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'P1 Start');
        plot(x2_initial(1), x2_initial(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'P2 Start');
        xlabel('X Position (m)');
        ylabel('Y Position (m)');
        title(sprintf('Particle Trajectories (e = %.2f)', e_for_dem));
        axis equal; % Ensure equal scaling for axes for correct perception of distances
        legend show;
        grid on; hold off;

        % Plot Particle Speeds over Time
        figure('Name', sprintf('Particle Speeds (e = %.2f)', e_for_dem));
        plot(t_vec_full, arrayfun(@norm, v1_history), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Particle 1 Speed'); hold on;
        plot(t_vec_full, arrayfun(@norm, v2_history), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Particle 2 Speed');
        xlabel('Time (s)');
        ylabel('Speed (m/s)');
        title(sprintf('Particle Speeds Over Time (e = %.2f)', e_for_dem));
        legend show;
        grid on; hold off;

        % Plot Overlap History (to visualize collision duration and penetration)
        figure('Name', sprintf('Overlap History (e = %.2f)', e_for_dem));
        plot(t_vec_full, overlap_history, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Overlap (\delta)');
        xlabel('Time (s)');
        ylabel('Overlap (m)');
        title(sprintf('Overlap Over Time (e = %.2f)', e_for_dem));
        grid on;

        % --- Extract Results for Problem 5 for return ---
        % These values will be passed back to the main function for display.
        
        % Initialize return values as NaN in case no overlap occurs
        v1_post_collision = NaN(1,3);
        v2_post_collision = NaN(1,3);
        duration_of_collision = NaN;
        contact_point_initial = NaN(1,3);
        contact_point_final = NaN(1,3);

        indices_overlap_positive = find(overlap_history > 0);
        if ~isempty(indices_overlap_positive)
            first_contact_idx = indices_overlap_positive(1);
            last_contact_idx = indices_overlap_positive(end);

            duration_of_collision = t_vec_full(last_contact_idx) - t_vec_full(first_contact_idx);
            
            % Find the index where particles separate (overlap becomes 0 after contact)
            idx_separation = find(overlap_history(last_contact_idx:end) == 0, 1, 'first') + last_contact_idx - 1;
            if isempty(idx_separation) || idx_separation >= num_steps
                % If separation not clearly detected by t_end, use velocities at the last step
                idx_separation = num_steps; 
            end
            
            v1_post_collision = v1_history(idx_separation, :);
            v2_post_collision = v2_history(idx_separation, :);

            % Initial contact point (conceptual for DEM, where overlap first becomes > 0)
            x1_at_first_contact = x1_history(first_contact_idx, :);
            x2_at_first_contact = x2_history(first_contact_idx, :);
            % Calls the nested function for Problem 3
            [~, n_hat_at_first_contact] = calculate_overlap_and_normal(x1_at_first_contact, d1, x2_at_first_contact, d2);
            contact_point_initial = x1_at_first_contact + r1 * n_hat_at_first_contact;

            % Contact point at separation (conceptual for DEM, where overlap returns to 0)
            x1_at_separation = x1_history(idx_separation, :);
            x2_at_separation = x2_history(idx_separation, :);
            % Calls the nested function for Problem 3
            [~, n_hat_at_separation] = calculate_overlap_and_normal(x1_at_separation, d1, x2_at_separation, d2);
            contact_point_final = x1_at_separation + r1 * n_hat_at_separation;

        else
            fprintf('  No overlap detected during simulation for e=%.2f. Check initial conditions or t_end.\n', e_for_dem);
        end
      end
        