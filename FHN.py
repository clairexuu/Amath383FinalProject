import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the FitzHugh-Nagumo model
def fhn_model(state, t, a, b, c):
    V, w = state
    dVdt = V * (a - V) * (V - 1) - w  # Assuming I = 0
    dwdt = b * V - c * w
    return [dVdt, dwdt]

# Parameters for three cases
cases = {
    "Case 1: Δ > 0 (Three Equilibria)": {"a": 0.1, "b": 0.01, "c": 0.02},
    "Case 2: Δ = 0 (Two Equilibria, Saddle-Node Bifurcation)": {"a": 0.25, "b": 0.01, "c": 0.02},
    "Case 3: Δ < 0 (One Equilibrium)": {"a": 0.5, "b": 0.01, "c": 0.02}
}

# Generate phase plane plots for each case
for case, params in cases.items():
    a, b, c = params["a"], params["b"], params["c"]

    # Create a grid for the phase plane
    V_vals = np.linspace(-0.5, 1.5, 20)
    w_vals = np.linspace(-0.2, 0.3, 20)
    V, W = np.meshgrid(V_vals, w_vals)

    # Compute the vector field
    dV = V * (a - V) * (V - 1) - W
    dW = b * V - c * W

    # Simulate trajectories
    t = np.linspace(0, 200, 2000)
    initial_conditions = [[-0.3, -0.1], [0.5, 0.1], [1.0, 0.2]]
    trajectories = [odeint(fhn_model, ic, t, args=(a, b, c)) for ic in initial_conditions]

    # Plot the phase plane
    plt.figure(figsize=(8, 6))
    plt.streamplot(V, W, dV, dW, color="gray", linewidth=0.7, density=1.2)

    # Plot nullclines
    V_range = np.linspace(-0.5, 1.5, 100)
    W_nullcline = V_range * (a - V_range) * (V_range - 1)
    V_nullcline = (b / c) * V_range

    plt.plot(V_range, W_nullcline, 'r', label="V-nullcline")
    plt.plot(V_range, V_nullcline, 'b', label="w-nullcline")

    # Solve for equilibrium points
    A = 1
    B = -(a + 1)
    C = a - (b / c)
    discriminant = B**2 - 4*A*C

    equilibria = [(0, 0)]  # The trivial equilibrium always exists
    if discriminant > 0:
        V1 = (-B + np.sqrt(discriminant)) / (2*A)
        V2 = (-B - np.sqrt(discriminant)) / (2*A)
        equilibria += [(V1, (b/c)*V1), (V2, (b/c)*V2)]
    elif discriminant == 0:
        V1 = -B / (2*A)
        equilibria.append((V1, (b/c)*V1))

    # Plot equilibrium points
    for eq in equilibria:
        plt.plot(eq[0], eq[1], 'go', markersize=8, label="Equilibrium" if eq == equilibria[0] else "")

    # Plot trajectories
    for traj in trajectories:
        plt.plot(traj[:, 0], traj[:, 1], 'k', alpha=0.8)

    plt.xlabel("Membrane Voltage (V)")
    plt.ylabel("Recovery Variable (w)")
    plt.title(f"Phase Plane of FitzHugh-Nagumo Model\n{case}")
    plt.legend()
    plt.grid()
    plt.show()