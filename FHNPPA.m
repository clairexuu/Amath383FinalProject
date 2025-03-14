clc; clear; close all;

% Parameters for the FHN model
I = 0; 
b = 0.01; 
a = 0.1; 
c = 0.1;

% Define the grid for the phase plane
v_range = linspace(-0.4, 1, 30);
w_range = linspace(-0.1, 0.3, 30);
[V, W] = meshgrid(v_range, w_range);

% FitzHugh-Nagumo equations
dV = V - (V.^3)/3 - W + I;  % Voltage equation
dW = b * (V + a - c * W);    % Recovery equation

% Plot phase portrait
figure;
quiver(V, W, dV, dW, 'k', 'AutoScale', 'on', 'AutoScaleFactor', 1);
hold on;

% Compute and plot nullclines
v_nullcline = V - (V.^3)/3 + I; % Set dV/dt = 0
w_nullcline = (V + a) / c;       % Set dW/dt = 0

plot(v_range, v_nullcline, 'k-', 'LineWidth', 1.5, 'DisplayName', 'V-nullcline');
plot(v_range, w_nullcline, 'k--', 'LineWidth', 1.5, 'DisplayName', 'W-nullcline');

% Find equilibrium points numerically by minimizing the squared error
equilibria = [];
for v_guess = linspace(-0.3, 0.9, 10)
    for w_guess = linspace(-0.1, 0.3, 10)
        eq_point = fminsearch(@(x) norm([x(1) - (x(1)^3)/3 - x(2) + I;
                                         b*(x(1) + a - c*x(2))]), [v_guess, w_guess]);
        if all(abs([eq_point(1) - (eq_point(1)^3)/3 - eq_point(2) + I;
                    b*(eq_point(1) + a - c*eq_point(2))]) < 1e-3) % Tolerance check
            equilibria = [equilibria; eq_point]; %#ok<AGROW>
        end
    end
end

% Remove duplicates
equilibria = unique(round(equilibria, 3), 'rows');

% Plot equilibrium points
plot(equilibria(:,1), equilibria(:,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

% Overlay a separatrix (approximate by numerically solving system)
tspan = [0, 100]; % Time span for integration
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
for i = 1:size(equilibria, 1)
    [~, sep_traj] = ode45(@(t, x) [x(1) - (x(1)^3)/3 - x(2) + I; 
                                   b * (x(1) + a - c * x(2))], ...
                          tspan, equilibria(i,:) + [0.1, 0.1], options);
    plot(sep_traj(:,1), sep_traj(:,2), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Separatrix');
end

% Annotate the plot
xlabel('Membrane voltage, V');
ylabel('Recovery, W');
title('FitzHugh-Nagumo Phase Portrait');
legend;
grid on;
hold off;
