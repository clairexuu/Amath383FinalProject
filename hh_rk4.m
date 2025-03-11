function [t,solution] = rk4_method(ode_RHS,t_range,num_steps,init_cond)
    % time steps
    step = (t_range(2) - t_range(1)) / num_steps;
    t = linspace(t_range(1), t_range(2), num_steps+1);
    solution = zeros(4, length(t));
    solution(:,1) = init_cond;

    % loop through each time step
    for j = 1:length(t)-1
        % calculate the four K values for Runge-Kutta method
        K1 = ode_RHS(t(j), solution(:,j));
        K2 = ode_RHS(t(j) + step/2, solution(:,j) + step*K1/2);
        K3 = ode_RHS(t(j) + step/2, solution(:,j) + step*K2/2);
        K4 = ode_RHS(t(j) + step, solution(:,j) + step*K3);

        % update step
        solution(:,j+1) = solution(:,j) + step/6 * (K1 + 2*K2 + 2*K3 + K4);
    end
end

% define the HH system
function dydt = ode_RHS(t, y)
    % unpack the state vectors
    V = y(1);  % membrane potential
    n = y(2);  % n gating variable
    m = y(3);  % m gating variable
    h = y(4);  % h gating variable

    % define constants
    g_K = 36;       % Potassium conductance (mS/cm^2)
    g_Na = 120;     % Sodium conductance (mS/cm^2)
    g_L = 0.3;      % Leak conductance (mS/cm^2)
    E_K = 12;      % Potassium equilibrium potential (mV)
    E_Na = -155;      % Sodium equilibrium potential (mV)
    E_L = -10.613;    % Leak equilibrium potential (mV)
    C_m = 1;        % Membrane capacitance (uF/cm^2)

    % alpha_x and beta_x functions for gating variables n, m, h
    alpha_n = 0.01 * (V + 10) / (exp((V + 10) / 10) - 1);
    beta_n = 0.125 * exp(V / 80);
    dndt = alpha_n * (1 - n) - beta_n * n;

    alpha_m = 0.1 * (V + 25) / (exp((V + 25) / 10) - 1);
    beta_m = 4 * exp(V / 18);
    dmdt = alpha_m * (1 - m) - beta_m * m;

    alpha_h = 0.07 * exp(V / 20);
    beta_h = 1 / (exp((V + 30) / 10) + 1);
    dhdt = alpha_h * (1 - h) - beta_h * h;

    % ion currents for V
    I_K = g_K * n^4 * (V - E_K);
    I_Na = g_Na * m^3 * h * (V - E_Na);
    I_L = g_L * (V - E_L);
    I = 10;  % (uA/cm^2)
    dVdt = (I - (I_K + I_Na + I_L)) / C_m;

    dydt = [dVdt; dndt; dmdt; dhdt];
end



%% Run runge_kutta_hh
% define time and num steps
t_range = [0 50]; % (ms)
num_steps = 2000; 
% when this is 1000 but euler's is 10000 there's a small difference

% set typical initial conditions and parameter I
V0 = -65;  % resting membrane potential (mV)
n0 = 0.3177;  % initial potassium gating variable
m0 = 0.0529;  % initial sodium activation gating variable
h0 = 0.5961;  % initial sodium inactivation gating variable
init_cond = [V0; n0; m0; h0];

% run runge-kutta-4
[t, solution] = rk4_method(@ode_RHS, t_range, num_steps, init_cond);



%% Plots
% plot change in V and gating variables
% V
figure;
subplot(2,1,1);
plot(t, solution(1,:), 'k');
xlabel('Time (ms)'); 
ylabel('Membrane Potential V (mV)');
title('Hodgkin-Huxley Model: Membrane Potential (RK4, numsteps=2000)');
grid on;

% gating variables
subplot(2,1,2);
plot(t, solution(2,:), 'r', t, solution(3,:), 'g', t, solution(4,:), 'b');
xlabel('Time (ms)'); 
ylabel('Gating Variables');
legend('n', 'm', 'h');
title('Gating Variables Over Time (RK4, numsteps=2000)');
grid on;

% plot gating variables vs. V
figure;
subplot(3,1,1);
plot(solution(2,:), solution(1,:), 'r');
xlabel('Potassium Activaton n');
ylabel('Memb Potential V (mV)'); 
title('n vs V (RK4, numsteps=2000)');
grid on;

subplot(3,1,2);
plot( solution(3,:), solution(1,:), 'g');
xlabel('Sodium Activation m');
ylabel('Memb Potential V (mV)'); 
title('m vs V (RK4, numsteps=2000)');
grid on;

subplot(3,1,3);
plot(solution(4,:), solution(1,:), 'b');
xlabel('Sodium Inactivation h');
ylabel('Memb Potential V (mV)'); 
title('h vs V (RK4, numsteps=2000)');
grid on;