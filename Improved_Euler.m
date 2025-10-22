clear all
clc
tstart = 0;
tend = 1;
dt = [10^-1 10^-2 10^-3 10^-4 10^-5 10^-6];
local_error = ones(1,length(dt));
global_error = ones(1,length(dt));
dt_count = length(dt)
f = @(y)-5*y;
for j = 1:dt_count;
    delta_t  = dt(j);
    t = [tstart:delta_t:tend];
    y = ones(1, length(t));
    y(1) = 1;
    for i = 1 : (length(t) - 1);
    y(i+1)= y(i) + delta_t*f(y(i));
    end
    figure(1)
    plot(t,y,'r');
    xlabel('time(t)');
    ylabel('height(y)');
    title('particle trajectory path')
    hold on
    y_analytic = y(1) .* exp(-5.*t);
    local_error(j) = abs( ( y(2) - y_analytic(2))/y_analytic(2)   );
    global_error(j) = abs( ( y(end) - y_analytic(end))/y_analytic(end)   );
end
    figure(1) 
    plot(t,y_analytic,'--');
    xlabel('t');
    ylabel('y');
    legend;
    hold on
    
figure(2)
hold on
loglog(dt,local_error,'r--','DisplayName','local, local_error');
hold on 
loglog(dt,global_error,'b--','DisplayName','global, global_error');
xlabel('dt');
ylabel('error');
legend;
function dydt = calculate_rhs_of_ODE(y)

  dydt = -5*y;
   
end