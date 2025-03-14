%% FitzHugh-Nagumo Model Simulation using Euler and RK4 Methods

clc; clear; close all;

%% Parameters
a = 0.7; 
b = 0.8; 
c = 0.08;
I = 0.5;  % External current

t_range = [0 50]; % Time range (ms)
num_steps = 2000; 
dt = (t_range(2) - t_range(1)) / num_steps;
t = linspace(t_range(1), t_range(2), num_steps+1);

% Initial conditions
V0 = -1;
w0 = 0;

%% Euler's Method
V_euler = zeros(1, num_steps+1);
w_euler = zeros(1, num_steps+1);
V_euler(1) = V0;
w_euler(1) = w0;

for j = 1:num_steps
    dVdt = V_euler(j) * (a - V_euler(j)) * (V_euler(j) - 1) - w_euler(j) + I;
    dwdt = b * V_euler(j) - c * w_euler(j);
    
    V_euler(j+1) = V_euler(j) + dt * dVdt;
    w_euler(j+1) = w_euler(j) + dt * dwdt;
end

%% RK4 Method
V_rk4 = zeros(1, num_steps+1);
w_rk4 = zeros(1, num_steps+1);
V_rk4(1) = V0;
w_rk4(1) = w0;

for j = 1:num_steps
    k1_V = dt * (V_rk4(j) * (a - V_rk4(j)) * (V_rk4(j) - 1) - w_rk4(j) + I);
    k1_w = dt * (b * V_rk4(j) - c * w_rk4(j));
    
    k2_V = dt * ((V_rk4(j) + k1_V/2) * (a - (V_rk4(j) + k1_V/2)) * ((V_rk4(j) + k1_V/2) - 1) - (w_rk4(j) + k1_w/2) + I);
    k2_w = dt * (b * (V_rk4(j) + k1_V/2) - c * (w_rk4(j) + k1_w/2));
    
    k3_V = dt * ((V_rk4(j) + k2_V/2) * (a - (V_rk4(j) + k2_V/2)) * ((V_rk4(j) + k2_V/2) - 1) - (w_rk4(j) + k2_w/2) + I);
    k3_w = dt * (b * (V_rk4(j) + k2_V/2) - c * (w_rk4(j) + k2_w/2));
    
    k4_V = dt * ((V_rk4(j) + k3_V) * (a - (V_rk4(j) + k3_V)) * ((V_rk4(j) + k3_V) - 1) - (w_rk4(j) + k3_w) + I);
    k4_w = dt * (b * (V_rk4(j) + k3_V) - c * (w_rk4(j) + k3_w));
    
    V_rk4(j+1) = V_rk4(j) + (k1_V + 2*k2_V + 2*k3_V + k4_V) / 6;
    w_rk4(j+1) = w_rk4(j) + (k1_w + 2*k2_w + 2*k3_w + k4_w) / 6;
end

%% Plot Results (Separate Figures)

% Euler's Method Plots
figure;
plot(t, V_euler, 'r', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Membrane Potential V');
title('Time vs Membrane Potential (Euler Method)');
grid on;

figure;
plot(t, w_euler, 'r', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Recovery Variable w');
title('Time vs Recovery Variable (Euler Method)');
grid on;

figure;
plot(V_euler, w_euler, 'r', 'LineWidth', 1.5);
xlabel('Membrane Potential V'); ylabel('Recovery Variable w');
title('Phase Portrait (Euler Method)');
grid on;

% RK4 Method Plots
figure;
plot(t, V_rk4, 'b', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Membrane Potential V');
title('Time vs Membrane Potential (RK4 Method)');
grid on;

figure;
plot(t, w_rk4, 'b', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Recovery Variable w');
title('Time vs Recovery Variable (RK4 Method)');
grid on;

figure;
plot(V_rk4, w_rk4, 'b', 'LineWidth', 1.5);
xlabel('Membrane Potential V'); ylabel('Recovery Variable w');
title('Phase Portrait (RK4 Method)');
grid on;
