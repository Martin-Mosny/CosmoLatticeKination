import numpy as np 
import matplotlib.pyplot as plt
import csv
from scipy.signal import find_peaks
from scipy.optimize import fsolve

# Read file
file_path = "/home/mosny/Desktop/cosmolattice/kination/average_scale_factor.txt"
data = np.loadtxt(file_path)

time = data[:, 0]
a_factor = data[:, 1]
a_factor_der = data[:, 2]

c = a_factor_der[0]
p_initial = 1

def equations(x, time, a_factor, a_factor_der):
    p, t0 = x

    equ1 = (1 + t_val / t0) ** p - a_factor

    equ2 = p * ( 1 + t_val / t0) ** (p - 1) / t0 - a_factor_der

    return [equ1, equ2]

p_solution = np.empty_like(time)
t0_solution = np.empty_like(time)

for i, (t_val, a_val, ad_val) in enumerate(zip(time, a_factor, a_factor_der)):

    initial_guess = [1.0, 1.0]
    sol = fsolve(equations, initial_guess, args=(t_val, a_val, ad_val))
    p_solution[i], t0_solution[i] = sol


# plt.plot(time, a_factor)
# plt.plot(time, (1+time/0.7)**0.33333333333)

plt.plot(time[8:], p_solution[8:])

plt.show()