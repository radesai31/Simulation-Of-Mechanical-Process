% Hand in 1 : Winnowing = Euler and RK4
% Name : Ridham Kishorbhai Bhuva
% Matr. No. : 245921
% Course : CEE

clear all;

% Dimensions and flow properties in SI unit
x0 = 0.50; 
y0 = 0.50;
xc = 0.55;
yc = -0.50;
h = 0.10;
u0 = 0.1;

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
N = 200;

% Time step
dt = 1/N;

% Time array and initial condition
t = zeros(1,N);
t(1) = 0;

% Initialization of (x,y) position, and intial conditions
x = zeros(N,2);
x(1,1) = x0; % Col 1 indicate Grain particle
x(1,2) = x0; % Col 2 indicate Chaff particle

y = zeros(N,2);
y(1,1) = y0; 
y(1,2) = y0;

% Initialization of (up_x,up_y) velocity and intial conditions 
up_x = zeros(N,2);
up_x(1,1) = 0;
up_x(1,2) = 0;

up_y = zeros(N,2);
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
        while y(i, n) > -0.5 % Bins are located at yc = -0.5

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

            % Functions of governing equations
            f_xx = @(t,up_x) up_x;
            f_xy = @(t,up_y) up_y;

            f_ux = @(t,up_x) F_dx/(rho_p*Vp); % Olny drag force in x direction
            f_uy = @(t,up_y) (1/(rho_p*Vp))*(F_g + F_b + F_dy);

            %% Euler Method calculation
            % Comment out one of the calculation for answer 1 and 2

            x(i+1,n) = x(i, n) + dt* f_xx(t(i),up_x(i,n));
            % particle is falling
            y(i+1,n) = y(i, n) - dt* f_xy(t(i),up_y(i,n)); 

            up_x(i+1,n) = up_x(i,n) + dt * f_ux(t(i),up_x(i,n));
            up_y(i+1,n) = up_y(i,n) + dt * f_uy(t(i),up_y(i,n));
        
            %% RK4 calculation
            
            kx1 = f_xx(t(i),up_x(i,n));
            kx2 = f_xx(t(i) + (0.5*dt),up_x(i,n) + (0.5*dt*kx1));
            kx3 = f_xx(t(i) + (0.5*dt),up_x(i,n) + (0.5*dt*kx2));
            kx4 = f_xx(t(i) + dt, up_x(i,n) + (dt*kx3));

            x(i+1,n) = x(i,n) + (1/6) * dt * (kx1 + 2*kx2 + 2*kx3 + kx4);

            ky1 = f_xy(t(i),up_y(i,n));
            ky2 = f_xy(t(i) + (0.5*dt),up_y(i,n) + (0.5*dt*ky1));
            ky3 = f_xy(t(i) + (0.5*dt),up_y(i,n) + (0.5*dt*ky2));
            ky4 = f_xy(t(i) + dt, up_y(i,n) + (dt*ky3));

            % particle is falling
            y(i+1,n) = y(i,n) - (1/6) * dt * (ky1 + 2*ky2 + 2*ky3 + ky4);
            

            kux1 = f_ux(t(i),up_x(i,n));
            kux2 = f_ux(t(i) + (0.5*dt), up_x(i,n) + (0.5*dt*kux1));
            kux3 = f_ux(t(i) + (0.5*dt), up_x(i,n) + (0.5*dt*kux2));
            kux4 = f_ux(t(i) + dt, up_x(i,n) + (dt*kux3));

            up_x(i+1,n) = up_x(i,n) + (1/6) * dt * (kux1 + 2*kux2 + 2*kux3 + kux4);

            kuy1 = f_uy(t(i),up_y(i,n));
            kuy2 = f_uy(t(i) + (0.5*dt), up_y(i,n) + (0.5*dt*kuy1));
            kuy3 = f_uy(t(i) + (0.5*dt), up_y(i,n) + (0.5*dt*kuy2));
            kuy4 = f_uy(t(i) + dt, up_y(i,n) + (dt*kuy3));

            up_y(i+1,n) = up_y(i,n) + (1/6) * dt * (kuy1 + 2*kuy2 + 2*kuy3 + kuy4);
         
            %% Step increase
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
