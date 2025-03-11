function [t, solution] = euler_method_hh(t_range, num_steps, init_cond, I)
    % define parameters
    g_K = 36;       % Potassium conductance (mS/cm^2)
    g_Na = 120;     % Sodium conductance (mS/cm^2)
    g_L = 0.3;      % Leak conductance (mS/cm^2)
    E_K = 12;      % Potassium equilibrium potential (mV)
    E_Na = -155;      % Sodium equilibrium potential (mV)
    E_L = -10.613;    % Leak equilibrium potential (mV)
    C_m = 1;        % Membrane capacitance (uF/cm^2)
    
    % time steps
    step = (t_range(2) - t_range(1)) / num_steps;
    t = linspace(t_range(1), t_range(2), num_steps+1);
    solution = zeros(4, num_steps+1);
    solution(:,1) = init_cond;
    
    % euler's method
    for j = 1:num_steps
        % grab previous iteration
        V = solution(1,j);
        n = solution(2,j);
        m = solution(3,j);
        h = solution(4,j);
        
        % compute alpha_x and beta_x
        alpha_n = 0.01 * (V + 10) / (exp((V + 10) / 10) - 1);
        beta_n = 0.125 * exp(V / 80);
        alpha_m = 0.1 * (V + 25) / (exp((V + 25) / 10) - 1);
        beta_m = 4 * exp(V / 18);
        alpha_h = 0.07 * exp(V / 20);
        beta_h = 1 / (exp((V + 30) / 10) + 1);
        
        % compute n, m, and h next derivatives
        dndt = alpha_n * (1 - n) - beta_n * n;
        dmdt = alpha_m * (1 - m) - beta_m * m;
        dhdt = alpha_h * (1 - h) - beta_h * h;
        
        % compute V next derivative
        I_K = g_K * n^4 * (V - E_K);
        I_Na = g_Na * m^3 * h * (V - E_Na); 
        I_L = g_L * (V - E_L);
        dVdt = (I - (I_K + I_Na + I_L)) / C_m;
        
        % euler update
        solution(1,j+1) = V + step * dVdt;
        solution(2,j+1) = n + step * dndt;
        solution(3,j+1) = m + step * dmdt;
        solution(4,j+1) = h + step * dhdt;
    end
end



%% Run euler_method_hh
% set time range and number of steps
t_range = [0 50]; % (ms)
num_steps = 2000;

% set typical initial conditions
V0 = -65;  % resting membrane potential (mV)
n0 = 0.3177;  % initial potassium gating variable
m0 = 0.0529;  % initial sodium activation gating variable
h0 = 0.5961;  % initial sodium inactivation gating variable
init_cond = [V0; n0; m0; h0];

% Set external current stimulus (constant or time-dependent)
I = 10;  % Applied current (uA/cm^2), change as needed

% run euler's method
[t, solution] = euler_method_hh(t_range, num_steps, init_cond, I);



%% Plots
% plot change in V and gating variables
% V
figure;
subplot(2,1,1);
plot(t, solution(1,:), 'k');
xlabel('Time (ms)'); 
ylabel('Membrane Potential V (mV)');
title('Hodgkin-Huxley Model: Membrane Potential (EULER, numsteps=2000)');
grid on;

% gating variables
subplot(2,1,2);
plot(t, solution(2,:), 'r', t, solution(3,:), 'g', t, solution(4,:), 'b');
xlabel('Time (ms)'); 
ylabel('Gating Variables');
legend('n', 'm', 'h');
title('Gating Variables Over Time (EULER, numsteps=2000)');
grid on;

% plot gating variables vs. V
figure;
subplot(3,1,1);
plot(solution(2,:), solution(1,:), 'r');
xlabel('Potassium Activaton n');
ylabel('Memb Potential V (mV)'); 
title('n vs V (EULER, numsteps=2000)');
grid on;

subplot(3,1,2);
plot( solution(3,:), solution(1,:), 'g');
xlabel('Sodium Activation m');
ylabel('Memb Potential V (mV)'); 
title('m vs V (EULER, numsteps=2000)');
grid on;

subplot(3,1,3);
plot(solution(4,:), solution(1,:), 'b');
xlabel('Sodium Inactivation h');
ylabel('Memb Potential V (mV)'); 
title('h vs V (EULER, numsteps=2000)');
grid on;