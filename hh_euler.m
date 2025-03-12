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