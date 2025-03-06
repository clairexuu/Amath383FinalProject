import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 0.1   # Controls shape of nullcline
b = 0.01  # Recovery variable kinetics
c = 0.02  # Recovery variable kinetics

# Define the quadratic equation coefficients
A = 1
B = -(a + 1)
C = a - (b / c)

# Compute the discriminant
D = B**2 - 4*A*C

# Compute equilibrium solutions
if D > 0:
    V1 = (-B + np.sqrt(D)) / (2*A)
    V2 = (-B - np.sqrt(D)) / (2*A)
    equilibria = [(0, 0), (V1, (b/c)*V1), (V2, (b/c)*V2)]
elif D == 0:
    V1 = -B / (2*A)
    equilibria = [(0, 0), (V1, (b/c)*V1)]
else:
    equilibria = [(0, 0)]  # Only the trivial equilibrium

# Print the equilibrium points
print("Equilibrium Solutions:")
for eq in equilibria:
    print(f"V* = {eq[0]:.4f}, w* = {eq[1]:.4f}")

# Plot equilibria
V_vals = np.linspace(-0.5, 1.5, 100)
W_nullcline = V_vals * (a - V_vals) * (V_vals - 1)
V_nullcline = (b / c) * V_vals

plt.figure(figsize=(8, 6))
plt.plot(V_vals, W_nullcline, 'r', label="V-nullcline")
plt.plot(V_vals, V_nullcline, 'b', label="w-nullcline")

for eq in equilibria:
    plt.plot(eq[0], eq[1], 'go', markersize=8)

plt.xlabel("Membrane Voltage (V)")
plt.ylabel("Recovery Variable (w)")
plt.title("Equilibrium Solutions of FitzHugh-Nagumo Model")
plt.legend()
plt.grid()
plt.show()
