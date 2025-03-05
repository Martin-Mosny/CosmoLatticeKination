import numpy as np 
import matplotlib.pyplot as plt
import csv
from scipy.signal import find_peaks
from scipy.optimize import fsolve
from scipy.misc import derivative

# Read file
file_path1 = "/home/mosny/Desktop/cosmolattice/kination/average_scale_factor.txt"
file_path2 = "/home/mosny/Desktop/cosmolattice/kination/average_energies.txt"
file_path3 = "/home/mosny/Desktop/cosmolattice/kination/average_scalar_0.txt"


data = np.loadtxt(file_path1)
data2 = np.loadtxt(file_path2)
data3 = np.loadtxt(file_path3)


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

def scipy_derivative(H: np.ndarray, t: np.ndarray):

    H_interp = lambda x: np.interp(x,t,H)
    dH_dt = np.array([derivative(H_interp, ti, dx = np.gradient(t).mean()) for ti in t])
    return dH_dt

p_solution = np.empty_like(time)
t0_solution = np.empty_like(time)

for i, (t_val, a_val, ad_val) in enumerate(zip(time, a_factor, a_factor_der)):

    initial_guess = [1.0, 1.0]
    sol = fsolve(equations, initial_guess, args=(t_val, a_val, ad_val))
    p_solution[i], t0_solution[i] = sol


# plt.plot(time, a_factor)
# plt.plot(time, (1+time/0.7)**0.33333333333)

# plt.plot(time[8:], p_solution[8:])

kinetic_energy = data2[:, 1]
gradient_energy = data2[:, 2]
potential_energy = data2[:,3]
total_energy = data2[:, 4]


def wave_average1(array, total_array):
    arrayUpper, _ = find_peaks(array, distance =  2)
    arrayLower, _ = find_peaks(-array, distance =  2)

    arrayPeaks = np.concatenate((arrayLower, arrayUpper))
    arrayPeaks.sort()

    index = 0
    average_array = np.array([])
    average_array_time = []
    value = 0
    for i in arrayPeaks:

        if index == 0:
            value = array[i] / total_array[i]
            index += 1
        elif index == 1:
            value += array[i] / total_array[i]
            average_array = np.append(average_array, (value / 2))
            average_array_time.append(i)
            index = 0
            value = 0
    return average_array, average_array_time

def wave_average2(array):
    arrayUpper, _ = find_peaks(array, distance =  2)
    arrayLower, _ = find_peaks(-array, distance =  2)

    arrayPeaks = np.concatenate((arrayLower, arrayUpper))
    arrayPeaks.sort()

    index = 0
    average_array = np.array([])
    average_array_time = []
    value = 0
    for i in arrayPeaks:

        if index == 0:
            value = array[i]
            index += 1
        elif index == 1:
            value += array[i]
            average_array = np.append(average_array, (value / 2))
            average_array_time.append(i)
            index = 0
            value = 0
    return average_array, average_array_time

average_kinetic, average_kinetic_time = wave_average1(kinetic_energy, total_energy)
average_gradient, average_gradient_time = wave_average1(gradient_energy, total_energy)

#plt.plot(time[8:], potential_energy[8:])

print(average_gradient_time)

#plt.plot(time, gradient_energy/total_energy)
#plt.plot(time, kinetic_energy/total_energy)
# plt.plot(time, potential_energy/total_energy)
# plt.plot(time[average_kinetic_time], average_kinetic)
# plt.plot(time[average_gradient_time], average_gradient)


scalar_vel = data3[:, 2]
scalar_vel_sqr = data3[:, 4]
scalar = data3[:, 1]
scalar_sqr =data3[:, 3]

H = data[:, 3]
dH_dt = scipy_derivative(H, time)

# plt.plot(time, scalar_vel * a_factor**2)

# Plotting H * a^2 is approx constant, showing it behaves like radiation
#plt.plot(time, H * a_factor**2)

x = scalar_vel / (4.2834 * np.sqrt(6) * H)
x2 = kinetic_energy / (3 * 4.2834**2 * H**2)
dx_dt = scipy_derivative(x, time)
dx2_dt = scipy_derivative(np.sqrt(x2),time)
# y = np.sqrt(0.05 * np.exp(- 0.8577 * scalar)) / (np.sqrt(3) * H)
y = np.sqrt(potential_energy) / (np.sqrt(3) * 4.2834 * H)
dy_dt = scipy_derivative(y, time)
z = np.sqrt(gradient_energy) /(np.sqrt(3) * 4.2834 * H)


#plt.plot(time, x2)
#plt.plot(time, x**2)
#plt.plot(time, y**2)
plt.plot(time, y**2)
plt.plot(time, x2 + y**2 + z**2)

x2_avg, x2_avg_time = wave_average2(x2)
z2_avg, z2_avg_time = wave_average2(z**2)
plt.plot(time[x2_avg_time], x2_avg)
plt.plot(time[z2_avg_time], z2_avg)

# Naive derivatives ignoring averaging difficulties
L = 0.8577
xprime = -3 * x + 4.2834 * L * np.sqrt(3/2) * y**2 + x * (1 + 2 * x2 - y**2)
yprime = -4.2834 * L * np.sqrt(3/2) * x * y + y * (1+2*x2 - y**2)
#plt.plot(time, H)
#plt.plot(time, dH_dt/ H**2)
#plt.plot(time, - (1 + 2* x**2 - y**2))

cross_term = 0.5 * np.sqrt(scalar_vel_sqr) * np.sqrt(scalar_sqr - scalar**2)
xcross = np.sqrt(scalar_vel_sqr - scalar_vel**2) /(np.sqrt(3) * 4.2834 * H)
picross = np.sqrt(scalar_vel_sqr - scalar_vel**2) / scalar_vel
phicross = np.sqrt(scalar_sqr - scalar**2)

# plt.plot(time, picross)
# plt.plot(time, x * y * phicross * np.sqrt(scalar_vel_sqr) /scalar_vel / 2 )

# plt.plot(time, dy_dt / H)
# plt.plot(time, yprime)
#plt.plot(time, dx_dt / H)
#plt.plot(time, dx2_dt / H)
#plt.plot(time, xprime)

#plt.plot(time, 1 + 2*x2 - y**2)

#plt.plot(time, potential_energy * a_factor**4)
#plt.plot(time, approx_potential_energy * a_factor**4)

yx = np.sqrt(0.05 * np.exp(-0.8577 * scalar)) / (np.sqrt(3)* 4.2834 * H)

# plt.plot(time, x2)
#plt.plot(time, 0.75*np.sqrt(6) * 4.2834 * np.log(a_factor))
# plt.plot(time, H)

plt.show()