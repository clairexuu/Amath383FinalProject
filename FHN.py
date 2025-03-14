import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# parameters
a, b, c = 0.1, 0.01, 0.1

# grid for the phase plane
V_vals = np.linspace(-0.5, 1.5, 20)
w_vals = np.linspace(-0.2, 0.3, 20)
V, W = np.meshgrid(V_vals, w_vals)

# vector field
dV = V * (a - V) * (V - 1) - W
dW = b * V - c * W

# nullclines
V_range = np.linspace(-0.5, 1.5, 100)
W_nullcline = V_range * (a - V_range) * (V_range - 1)
V_nullcline = (b / c) * V_range

# equilibrium points
V_sym = sp.Symbol('V')
eq_expr = V_sym * (a - V_sym) * (V_sym - 1) - (b/c) * V_sym

V_solutions = sp.solve(eq_expr, V_sym)  # Solve for V
equilibria = [(float(V), float((b/c) * V)) for V in V_solutions]  # Compute (V*, w*)

# Plot the phase plane
plt.figure(figsize=(8, 6))
plt.streamplot(V, W, dV, dW, color="gray", linewidth=0.7, density=1.2)

# Plot nullclines
plt.plot(V_range, W_nullcline, 'r', label="V-nullcline")
plt.plot(V_range, V_nullcline, 'b', label="w-nullcline")

# Plot equilibrium points
for eq in equilibria:
    plt.plot(eq[0], eq[1], 'go', markersize=8)

plt.xlabel("Membrane Voltage (V)")
plt.ylabel("Recovery Variable (w)")
plt.title(f"Phase Plane of FitzHugh-Nagumo Model\n(a={a}, b={b}, c={c})")
plt.legend()
plt.grid()
plt.show()

# Print the equilibrium points
print("Equilibrium points:", equilibria)