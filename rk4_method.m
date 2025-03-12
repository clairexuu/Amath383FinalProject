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