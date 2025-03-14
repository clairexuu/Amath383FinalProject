import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.optimize import fsolve

# Define Hodgkin-Huxley parameters
g_K = 36  # Potassium conductance (mS/cm^2)
g_Na = 120  # Sodium conductance (mS/cm^2)
g_L = 0.3  # Leak conductance (mS/cm^2)
E_K = 12  # Potassium equilibrium potential (mV)
E_Na = -155  # Sodium equilibrium potential (mV)
E_L = -10.613  # Leak equilibrium potential (mV)
C_m = 1  # Membrane capacitance (uF/cm^2)

# Define symbolic variables
V_sym, n_sym, m_sym, h_sym = sp.symbols('V n m h')

# Steady-state activation/inactivation functions
α_n = 0.01 * (10 - V_sym) / (sp.exp((10 - V_sym) / 10) - 1)
β_n = 0.125 * sp.exp(-V_sym / 80)

α_m = 0.1 * (25 - V_sym) / (sp.exp((25 - V_sym) / 10) - 1)
β_m = 4 * sp.exp(-V_sym / 18)

α_h = 0.07 * sp.exp(-V_sym / 20)
β_h = 1 / (sp.exp((30 - V_sym) / 10) + 1)

n_inf = α_n / (α_n + β_n)
m_inf = α_m / (α_m + β_m)
h_inf = α_h / (α_h + β_h)

# Convert symbolic expressions to numerical functions using lambdify
hh_n_inf = sp.lambdify(V_sym, n_inf, 'numpy')
hh_m_inf = sp.lambdify(V_sym, m_inf, 'numpy')
hh_h_inf = sp.lambdify(V_sym, h_inf, 'numpy')


# Define numerical function for equilibrium computation
def hh_equilibrium(V):
    """Compute equilibrium by setting dV/dt = 0."""
    n_eq = hh_n_inf(V)
    m_eq = hh_m_inf(V)
    h_eq = hh_h_inf(V)

    I_K = g_K * (n_eq ** 4) * (V - E_K)
    I_Na = g_Na * (m_eq ** 3) * h_eq * (V - E_Na)
    I_L = g_L * (V - E_L)

    return I_K + I_Na + I_L  # Set to zero for equilibrium condition


# Find equilibrium points numerically using fsolve
V_initial_guesses = np.linspace(-100, 50, 5)  # Different initial guesses
V_equil_vals = []
for guess in V_initial_guesses:
    V_eq = fsolve(hh_equilibrium, guess)[0]
    if np.isreal(V_eq) and V_eq not in V_equil_vals:
        V_equil_vals.append(V_eq)

# Compute corresponding (n, m, h) values at equilibrium points
equilibrium_points = []
for v_eq in V_equil_vals:
    eq_n = hh_n_inf(v_eq)
    eq_m = hh_m_inf(v_eq)
    eq_h = hh_h_inf(v_eq)
    equilibrium_points.append((v_eq, eq_n, eq_m, eq_h))

# Define HH current equations
I_K = g_K * (n_sym ** 4) * (V_sym - E_K)
I_Na = g_Na * (m_sym ** 3) * h_sym * (V_sym - E_Na)
I_L = g_L * (V_sym - E_L)

# Hodgkin-Huxley membrane equation
dV_dt = (-(I_K + I_Na + I_L)) / C_m

# Gating variable equations
dn_dt = α_n * (1 - n_sym) - β_n * n_sym
dm_dt = α_m * (1 - m_sym) - β_m * m_sym
dh_dt = α_h * (1 - h_sym) - β_h * h_sym

# Define state vector and compute Jacobian matrix
state_vars = [V_sym, n_sym, m_sym, h_sym]
HH_system = [dV_dt, dn_dt, dm_dt, dh_dt]
J = sp.Matrix(HH_system).jacobian(state_vars)

# Compute eigenvalues at equilibrium points
stability_results = {}
for eq_point in equilibrium_points:
    J_eval = J.subs({V_sym: eq_point[0], n_sym: eq_point[1], m_sym: eq_point[2], h_sym: eq_point[3]})
    eigenvalues = J_eval.eigenvals()
    stability_results[eq_point] = [ev.evalf() for ev in eigenvalues]

# Print equilibrium points and stability analysis
print("\nEquilibrium Points and Stability Analysis:")
for eq, eig_vals in stability_results.items():
    print(f"\nEquilibrium Point: (V, n, m, h) = {eq}")
    print("Eigenvalues:", eig_vals)

# Plot phase plane: Steady-state n_inf(V) and equilibrium points
V_vals = np.linspace(-100, 50, 100)
n_vals = [float(n_inf.subs(V_sym, v)) for v in V_vals]

plt.figure(figsize=(8, 6))
plt.plot(V_vals, n_vals, 'r', label=r'$n_{\infty}(V)$')
for v_eq in V_equil_vals:
    plt.axvline(x=v_eq, color='g', linestyle='--', label=f'Equilibrium: V={v_eq:.2f} mV')

plt.xlabel("Membrane Voltage (V)")
plt.ylabel("Steady-State Activation Variable n")
plt.title("Phase Plane Analysis of Hodgkin-Huxley Model")
plt.legend()
plt.grid()
plt.show()