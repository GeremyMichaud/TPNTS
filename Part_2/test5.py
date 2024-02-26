from scipy.stats import rv_continuous, maxwell
import numpy as np
import matplotlib.pyplot as plt
from vpython import *

mass_electron = 9E-31
k = 1.4E-23  # Boltzmann constant # TODO: changer pour une constante de Boltzmann en eV/K :)
t = 300  # around room temperature
class mb_speed(rv_continuous):
    def _pdf(self, x):
        return x * mass_electron / (k*t) * np.exp(-(mass_electron * x ** 2) / (2 * k * t))

distribution = mb_speed(a=0)
test = distribution.rvs(size=10000)

plt.hist(test, bins=100)
plt.show()