from numpy import *
from math import *
import itertools
import csv
from scipy import *
import scipy.integrate as scintegrate
from scipy.integrate import quad 
import scipy.special as sc
import sdeint as sd


Num_part = 10000
gamma = -0.4
D_x = 0.1
D_y = 0.1
L = 1000.0
a = 5
T = 50000

#if YES the harmonics are sinus, else they are cosine
SIN = "NO"

(tspan, step) = linspace(0.0, T, 25000, retstep=True)
N_max = int(ceil(L/(2*a)))
print(N_max)

x0 = 0.0
Y_init_cdts= linspace(-L/2,L/2,Num_part)
#Y_init_cdts= [-0.1,123]


#phi = zeros(N_max+1)
phi = random.rand(N_max+1)*2*math.pi
if SIN == "YES":
    write_file_phases = 'numerics/phases_L_' + str(L) + '_a_' + str(a) + '_gamma_' + str(gamma) + '_D_' + str(D_x) + '_sin_SR2_' + str(Num_part) + 'part_1.txt'
else:
    write_file_phases = 'numerics/phases_L_' + str(L) + '_a_' + str(a) + '_gamma_' + str(gamma) + '_D_' + str(D_x) + '_cos_SR2_' + str(Num_part) + 'part_1.txt'
file_phases = open(write_file_phases,'w')
for n in range(0,len(phi)):
    s = str(n) + '  '  + str(phi[n])  +'\n'
    file_phases.write(s)
file_phases.close()





def drift_Y(y, t):
    return 0

def noise_Y_ampl(y, t):
    return sqrt(2*D_y)

def noise_X_ampl(x, t):
    return sqrt(2*D_x)

def vel_field(x, t):
    if trajectory_Y[int(t/step)] > L/2:
        trajectory_Y[int(t/step)] =  trajectory_Y[int(t/step)] - L
    if trajectory_Y[int(t/step)] < -L/2:
        trajectory_Y[int(t/step)] =  trajectory_Y[int(t/step)] + L
    sum_k = 0.0
    for n in range(1,N_max+1):
        kn = 2*math.pi*n/L
        U_kn = power(kn,gamma/2)
        if SIN == "YES":
            sum_k  =  sum_k + U_kn*sin(kn*trajectory_Y[int(t/step)]+phi[n])
        else:
            sum_k  =  sum_k + U_kn*cos(kn*trajectory_Y[int(t/step)]+phi[n])
    return sum_k


if SIN == "YES":
    write_file_trajectory = 'numerics/trajectories_L_' + str(L) + '_a_' + str(a) + '_gamma_' + str(gamma) + '_D_' + str(D_x) + '_sin_SR2_' + str(Num_part) + 'part_1.csv'
else:
    write_file_trajectory = 'numerics/trajectories_L_' + str(L) + '_a_' + str(a) + '_gamma_' + str(gamma) + '_D_' + str(D_x) + '_cos_SR2_' + str(Num_part) + 'part_1.csv'
with open(write_file_trajectory, 'w', newline='') as file:
        fieldnames = ['time', 'X_position', 'Y_position', 'Y0']
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for y0 in Y_init_cdts:
            #trajectory_Y = sd.itoint(drift_Y, noise_Y_ampl, y0, tspan)
            #trajectory_Y = sd.stratint(drift_Y, noise_Y_ampl, y0, tspan)
            #trajectory_Y = sd.itoEuler(drift_Y, noise_Y_ampl, y0, tspan)
            #trajectory_Y = sd.stratHeun(drift_Y, noise_Y_ampl, y0, tspan)
            trajectory_Y = sd.itoSRI2(drift_Y, noise_Y_ampl, y0, tspan)
            #trajectory_Y = sd.stratSRS2(drift_Y, noise_Y_ampl, y0, tspan)
            #trajectory_Y = sd.stratKP2iS(drift_Y, noise_Y_ampl, y0, tspan)
            
            #trajectory_X = sd.itoint(vel_field, noise_X_ampl, x0, tspan)
            #trajectory_X = sd.stratint(vel_field, noise_X_ampl, x0, tspan)
            #trajectory_X = sd.itoEuler(vel_field, noise_X_ampl, x0, tspan)
            #trajectory_X = sd.stratHeun(vel_field, noise_X_ampl, x0, tspan)
            trajectory_X = sd.itoSRI2(vel_field, noise_X_ampl, x0, tspan)
            #trajectory_X = sd.stratSRS2(vel_field, noise_X_ampl, x0, tspan)
            #trajectory_X = sd.stratKP2iS(vel_field, noise_X_ampl, x0, tspan)
            #print(trajectory_X)

            for t in tspan:
                writer.writerow({'time': t, 'X_position': trajectory_X[int(t/step)][0],'Y_position': trajectory_Y[int(t/step)][0], 'Y0': y0})
