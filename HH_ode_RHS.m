% define the HH system
function dydt = HH_ode_RHS(t, y)
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

