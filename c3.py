import numpy as np

import matplotlib.pyplot as plt

# Constants
c = 3e8        # Speed of light in m/s
n = 1.5        # Refractive index (e.g., glass)
L = 2.0        # Distance between two points on x-axis (meters)
y = 1.0        # Height above x-axis (meters)

# x values from 0 to L
x = np.linspace(0, L, 500)

# Fermat's principle equation for time t
t = (np.sqrt(x**2 + y**2) / (c / n)) + (np.sqrt((L - x)**2 + y**2) / (c / n))

# find minimum
i = np.argmin(t)
x_min = x[i]
t_min = t[i]

# because of linspace, there prolly isnt an x = 0, but more likely x = 0.999, x = 1.001 etc meaning that the minimum is not found at x = 1
# therefore need rounding
print(np.round(x_min))
print(L/2)

if np.round(x_min) == L/2:
    print("Minimum occurs at the midpoint of the path, confirming symmetry.")

# Plot
plt.figure(figsize=(8, 5))
plt.plot(x, t, label='t(x) by Fermat\'s Principle')

plt.plot(x_min, t_min, marker='o')  # Mark the minimum point

plt.xlabel('x (m)')
plt.ylabel('t (s)')
plt.title('Parabola: x vs t (Fermat\'s Principle)')
plt.legend()
plt.grid(True)
plt.show()