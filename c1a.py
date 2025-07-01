import math
import matplotlib.pyplot as plt

'''
Brief: 
Use the seillmeier formula to plot the refractive index (n) of crown glass vs wavelength (400-800nm)
'''

c = 300000000 # speed of light
n = 1.5
cs = c/n

'''
sellmeier equation with 3 terms:
n^2 = 1 + (b1*lambda2)/(lambda2-c1) + (b2*lambda2)/(lambda2-c2) (b3*lambda2)/(lambda2-c3) 
'''

# sellmeier coeffs from wikipedia 
b = [1.03961212, 0.231792344, 1.01046945]
c = [0.00600069867, 0.0200179144, 103.560653] # milimetres squared

wavelengths = range(400, 801)
n = [math.sqrt(1 + (b[0]*(wavelength**2))/((wavelength**2)-c[0]) + (b[1]*(wavelength**2))/((wavelength**2)-c[1]) + (b[2]*(wavelength**2))/((wavelength**2)-c[2])) for wavelength in [w/1000 for w in wavelengths]] # do w/1000 as sellmeier expects micrometers not nanometres

fig, ax = plt.subplots()
ax.plot(wavelengths, n)
plt.show()