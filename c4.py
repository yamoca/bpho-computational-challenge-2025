import numpy as np
from matplotlib import pyplot as plt


# Constants
c = 3e8        # Speed of light in m/s
n1, n2 = 1.0, 1.5  # Refractive indices
c1, c2 = c/n1, c/n2  # Speed of light in media
L = 2.0        # Distance between two points on x-axis (meters)
y = 1.0        # Height above x-axis (meters)

# x values from 0 to L
x = np.linspace(0, L, 500)

# Fermat's principle equation for time t
t = (np.sqrt(x**2 + y**2) / c1) + (np.sqrt((L - x)**2 + y**2) / c2)
# previous equation changes c / n to c1 because speed of light in medium is c / n which is c1

# find minimum
i = np.argmin(t)
x_min = x[i]
t_min = t[i]

# initially i tried to::
# find theta and phi at minimum
# as x = ytan theta and l-x = ytan phi: (remember that t is our y axis)
theta = np.arctan(y / x_min)
phi = np.arctan(y / (L - x_min))
print(theta, phi)

# but really i wanted to find show that sin theta / c1 = sin phi / c2
# therefore i needed sin not tan: and instead of finding theta and then doing sin theta, i can just use sin directly to increase accuracy
# if we look at diagram, we see that as sin = opposite / hypotenuse, sin theta = x / sqrt(x^2 + y^2) and sin phi = (L - x) / sqrt((L - x)^2 + y^2)
sin_theta = x_min / np.sqrt(x_min**2 + y**2)
sin_phi = (L - x_min) / np.sqrt((L - x_min)**2 + y**2)
# now we can show that sin theta / c1 = sin phi / c2
if np.round(sin_theta / c1) == np.round(sin_phi / c2): # again floating point accuracy so convert
    print("Fermat's principle holds: sin(theta)/c1 = sin(phi)/c2")

print(sin_theta/c1, sin_phi/c2)
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