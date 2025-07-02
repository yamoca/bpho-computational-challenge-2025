import math
from matplotlib import pyplot as plt
frequencyaxis = range(405, 790)
frequency = [f*(10**12) for f in frequencyaxis]
# n = [math.sqrt(1+math.sqrt(1/(1.731-0.261*(f/10**15)**2))) for f in frequency]
n = [math.sqrt(1+math.sqrt(1/(1.731-0.261*(f/10**15)**2))) for f in frequency]

fig, ax = plt.subplots()
ax.plot(frequencyaxis, n)
plt.xlabel("frequency THz")
plt.ylabel("refractive index")
plt.show()


# initla problems thinking icould model the whole range from 0 to 1500 thz, missing the key information that the formula was only valid for visible light spectrum
# after fixing this, slightly bamboozled by converting between hz and thz and back again, so i seperated out axis numbers and numbers for use 