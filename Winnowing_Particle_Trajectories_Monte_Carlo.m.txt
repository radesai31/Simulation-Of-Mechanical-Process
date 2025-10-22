% Hand in 2 : Winnowing = monte Carlo
% Name : Ridham Kishorbhai Bhuva
% Matr. No. : 245921
% Course : CEE

%% Task 1 and 2
% please copy separate task for task 3 in new file and run it

clear all;

% Dimensions and flow properties in SI unit
x0 = 0.50; 
y0 = 0.50;
xc = 0.55;
yc = -0.50;
h = 0.10;
u0 = 0.2;

% Fluid properties (air) in SI unit
rho_f = 1.2;
mu_f = 1.8e-5;

% Grain particle properties in SI unit
rho_grain = 750;
dp_grain = 2.5e-3;

% Chaff particle properties in SI unit
rho_chaff = 50;
dp_chaff = 3.25e-3;

% Gravity constant
g = 9.81;

% Number of steps
M = 100;

% Number of samples
N = 10;

% Time step
dt = 0.5/M;

% Time array and initial condition
t = zeros(1,M+1);
t(1) = 0;

% Initialization of (x,y) position, and intial conditions
x = zeros(M+1,2);
x(1,1) = x0; % Col 1 indicate Grain particle
x(1,2) = x0; % Col 2 indicate Chaff particle

y = zeros(M+1,2);
y(1,1) = y0; 
y(1,2) = y0;

% Initialization of (up_x,up_y) velocity and intial conditions 
up_x = zeros(M+1,2);
up_x(1,1) = 0;
up_x(1,2) = 0;

up_y = zeros(M+1,2);
up_y(1,1) = 0;
up_y(1,2) = 0;

% Function of fluid velocity
f_uf = @(x,y) 6.2 * u0 * sqrt(h/x) * exp((-50)*((y^2)/(x^2)));

% Function of vector magnitude
f_mag = @(x,y) sqrt((x^2)+(y^2));

for n = 1:2 % n=1 for Grain particle and n=2 for Chaff particle
    if n == 1
        dp = dp_grain;
        rho_p = rho_grain;
    else
        dp = dp_chaff;
        rho_p = rho_chaff;
    end
    
    % volume of particles
    Vp = (pi*(dp^3))/6;
    
    i = 1;
    while y(i, n) > yc % Bins are located at yc = -0.5

          uf = f_uf(x(i,n),y(i,n));

            % Magnitude of uf minus up ||uf-up||
            % uf and up are perpendicular, uf vecotr only in x direction 
            uf_minus_up = f_mag((uf - up_x(i,n)),(0 - up_y(i,n)));

            % Re of particle
            Re_p = (rho_f * uf_minus_up * dp)/ mu_f;

        % Drag coefficinet Cd
        if Re_p < 800
            Cd = (24/Re_p) * ( 1 + 0.15*(Re_p^0.687) );
        else
            Cd = 0.44;
        end

            % Gravitational force
            F_g = rho_p * Vp * g;

            % Buoyancy force
            F_b = -(rho_f * Vp * g);

            % Drag force
            % In x direction
            F_dx = 0.5*pi*dp*dp*rho_f*Cd* uf_minus_up * (uf - up_x(i,n));
            
            % In y direction
            F_dy = 0.5*pi*dp*dp*rho_f*Cd* uf_minus_up * (0 - up_y(i,n));

            % Euler and monte carlo combine calculation for treajectory
            MCsum = 0;

            for j = 1:N
                trand = t(i) + rand() * dt;
                MCsum = MCsum + trand;
            end
            x(i+1,n) = x(i,n) + dt * up_x(i,n) * MCsum / N;
            y(i+1,n) = y(i,n) - dt * up_y(i,n) * MCsum / N;

            up_x(i+1,n) = up_x(i,n) + dt * (F_dx) * MCsum / (rho_p * Vp * N);
            up_y(i+1,n) = up_y(i,n) + dt * (F_g + F_b + F_dy) * MCsum / (rho_p * Vp * N);

            % Step increase
            t(i+1) = t(i) + dt;
            i = i+1;
    end
end
        
x_grain = x(x(:,1) ~= 0,1); %without zeros
y_grain = y(x(:,1) ~= 0,1);

x_chaff = x(:,2);
y_chaff = y(:,2);

figure(1)
plot(x_grain, y_grain, 'b-', 'LineWidth',2);
hold on
plot(x_chaff, y_chaff, 'r-', 'LineWidth',2);
hold off
xlabel('x (Bins)');
ylabel('y')
xline(xc, 'k--')
legend('Grain particle', 'Chaff particle');
title('Trajectory of particle');

% task 1 and 2 END!!!

%% TASK 3
% Please copy separate section in new file to run the task

clear all;

% Dimensions and flow properties in SI unit
x0 = 0.50; 
y0 = 0.50;
xc = 0.55; % X-coordinate of bin separation
yc = -0.50;
h = 0.10;
u0 = 0.2;

% Fluid properties (air) in SI unit
rho_f = 1.2;
mu_f = 1.8e-5;

% Gravity constant
g = 9.81;

% Number of steps
M = 100;

% 1000 particle samples
N = 1000;

% sample for monte carlo
N2 = 10;

% Time step
dt = 0.5/M;

% Time array and initial condition
t = zeros(1,M+1);
t(1) = 0;

% Initialization of (x,y) position, and intial conditions
x = zeros(M+1,2);
x(1,1) = x0; % Col 1 indicate Grain particle
x(1,2) = x0; % Col 2 indicate Chaff particle

y = zeros(M+1,2);
y(1,1) = y0; 
y(1,2) = y0;

% Initialization of (up_x,up_y) velocity and intial conditions 
up_x = zeros(M+1,2);
up_y = zeros(M+1,2);

% Angle between up and x-axis in degree
a1 = -95;
b1 = -85;
theta = a1 + (b1 - a1).*rand(1,N);

% Uniformaly distributed light particle diameter
a = 0.002;
b = 0.005;
dp_chaff = a + (b-a).*rand(1,N);

% Normally distributed heavy particle diameter
sigma_d = 0.001;
dp_mean = 0.0025;
dp_grain = randn(1,N).* sigma_d + dp_mean;

% Density of heavy particle
rho_grain = 750;

% Normally distributed density of light particle
sigma_rho = 20;
rho_mean = 50;
rho_chaff = randn(1,N).*sigma_rho + rho_mean;

% Initial zero particle in both bins
bin1 = 0;
bin2 = 0;

% Selecting Grain or Chaff particle 
% Plese use n = 1 for grain and n = 2 for Chaff
n = 1;

% Function of fluid velocity
f_uf = @(x,y) 6.2 * u0 * sqrt(h/x) * exp((-50)*((y^2)/(x^2)));

% Function of vector magnitude
f_mag = @(x,y) sqrt((x^2)+(y^2));

% Terminal velocity function
function ut = ut(rho_p, d_p, rho_f, g)
    ut = sqrt((4/3) * (rho_p - rho_f) .* g .* d_p ./ (3 * rho_f));
end

for n = 1:2 % 1000 grain particle and 1000 chaff particles
    for m = 1:N
        
            % dp and rho as per grain and chaff particle
                if n == 1
                    dp = dp_grain(m); % for n = 1
                    rho_p = rho_grain; % for n = 1
                else
                    dp = dp_chaff(m); % for n = 2
                    rho_p = rho_chaff(m); % for n = 2
                end
                    
            % Initial velocity as per terminal velocity and angle between
            % x-axis and velocity vector
            up_x(1,1) = abs(ut(rho_p,dp,rho_f,g) * cosd(180 - theta(m)));
            up_y(1,1) = abs(ut(rho_p,dp,rho_f,g) * cosd(theta(m) + 90));

            up_x(1,2) = abs(ut(rho_p,dp,rho_f,g) * cosd(180 - theta(m)));
            up_y(1,2) = abs(ut(rho_p,dp,rho_f,g) * cosd(theta(m) + 90));
            
            % volume of particles
            Vp = (pi*(dp^3))/6;
    
            i = 1;
            while y(i, n) > yc % Bins are located at yc = -0.5
    
              uf = f_uf(x(i,n),y(i,n));
    
                % Magnitude of uf minus up ||uf-up||
                % uf and up are perpendicular, uf vecotr only in x direction 
                uf_minus_up = f_mag((uf - up_x(i,n)),(0 - up_y(i,n)));
    
                % Re of particle
                Re_p = (rho_f * uf_minus_up * dp)/ mu_f;
    
            % Drag coefficinet Cd
            if Re_p < 800
                Cd = (24/Re_p) * ( 1 + 0.15*(Re_p^0.687) );
            else
                Cd = 0.44;
            end

                % Gravitational force
                F_g = rho_p * Vp * g;

                % Buoyancy force
                F_b = -(rho_f * Vp * g);

                % Drag force
                % In x direction
                F_dx = 0.5*pi*dp*dp*rho_f*Cd* uf_minus_up * (uf - up_x(i,n));
            
                % In y direction
                F_dy = 0.5*pi*dp*dp*rho_f*Cd* uf_minus_up * (0 - up_y(i,n));

                % Euler and monte carlo combine calculation for treajectory
                MCsum = 0;

                for k = 1:N2
                    trand = t(i) + rand() * dt;
                    MCsum = MCsum + trand;
                end
                    x(i+1,n) = x(i,n) + dt * up_x(i,n) * MCsum / N2;
                    y(i+1,n) = y(i,n) - dt * up_y(i,n) * MCsum / N2;

                    up_x(i+1,n) = up_x(i,n) + dt * (F_dx) * MCsum / (rho_p * Vp * N2);
                    up_y(i+1,n) = up_y(i,n) + dt * (F_g + F_b + F_dy) * MCsum / (rho_p * Vp * N2);

            % Step increase
            t(i+1) = t(i) + dt;
            i = i+1;

            end
        
                % Particle final position in bin1 or bin2
                if x(i,n) < xc
                    bin1 = bin1 + 1;
                end
            
                if x(i,n) > xc
                    bin2 = bin2 + 1;
                end
    end
end

Proportion_of_particles_in_Bin_1 = bin1/20 % in %
Proportion_of_particles_in_Bin_2 = bin2/20 % in %


% task 3 END!!!
