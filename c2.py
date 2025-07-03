import numpy as np
from matplotlib import pyplot as plt

u = np.arange(20, 60, 5)
v = np.array([65.5, 40.0, 31.0, 27.0, 25.0, 23.1, 21.5, 20.5])

# u_recip = [1/u for u in u]
# v_recip = [1/v for v in v]
# dont have to do this in numpy jut do this
u_recip = 1/u
v_recip = 1/v

best_fit_line_coeffs = np.polyfit(u_recip, v_recip, 1)
best_fit_line = np.poly1d(best_fit_line_coeffs)
x = np.linspace(min(u_recip), max(u_recip), 100)
line_equation = f"y = {best_fit_line_coeffs[0]}x + {best_fit_line_coeffs[1]}" 

fig, ax = plt.subplots()
ax.plot(u_recip, v_recip, ".")
ax.plot(x, best_fit_line(x), "-", label=line_equation)
ax.legend()
plt.xlabel("1/u (cm)")
plt.ylabel("1/v (cm)")
plt.show()

# points = np.array([x_data, y_data]).T

# learning about numpy e.g numpy arrays vs python lists, speed stuff. .T (transpose) very useful
# arrange is like range but with numpy array

# do some fancy stats stuff to assess "veracity" e.g r^2 stuff or pmcc and strong pos correlation or whatever