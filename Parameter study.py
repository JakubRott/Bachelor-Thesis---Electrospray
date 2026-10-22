import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import qmc
import pandas as pd


#Definition of physical parameters of deionaze water
rho_water = 1000 # Density of water in kg/m^3
gamma = 0.072 # Surface tension of water in N/m
K = 3.6e-4 # Electrical conductivity of water in S/m
epsilon_r = 80.4 # Relative permittivity of water
nu = 1e-3 # Dynamic viscosity of water in Pa.s
epsilon_0 = 8.854e-12 # Vacuum permittivity in F/m

transfer_ml_h_to_m3_s = 1e-6/3600 # Conversion factor from ml/h to m^3/s

#bound_parameters
boundary_parameters = {
    'd_in': (20e-6, 250e-6),                                        #in m however its defined in the paper as 20-250 um
    'd_out': (50e-6, 500e-6),                                       #in m however its defined in the paper as 50-500 um
    'R_corner': (10e-6, 50e-6),                                     #in m however its defined in the paper as 10-50 um
    'H': (3e-3,35e-3),                                              #in m however its defined in the paper as 3-35 mm
    'Q': (1e-3*transfer_ml_h_to_m3_s, 5*transfer_ml_h_to_m3_s),     #in m^3/s however its defined in the paper as 1-5 ml/h
    'V': (1.5e3,15e3)                                               #V however its defined in the paper as 1.5-15 kV 
}

resolution = {
    'd_in':     1e-6,                           # 1 µm   (lithography/etching tolerance)
    'd_out':    2e-6,                           # 2 µm
    'R_corner': 2e-6,                           # 2 µm
    'H':        0.5e-3,                         # 0.5 mm (positioning stage)
    'Q':        0.01 * transfer_ml_h_to_m3_s,   # 0.01 ml/h (syringe pump resolution)
    'V':        50,                             # 50 V   (HV supply setting)
}


Q0 = (epsilon_0*gamma)/(rho_water*K) #Characteristic flow rate in m^3/s
R0 = ((gamma*epsilon_0**2)/(rho_water*K**2))**(1/3) #Characteristic radius in m

def Rj(Q):
    return epsilon_r**(1/6)*R0*(Q/Q0)**(1/3)

def V_onset(d_out, H):
    return np.sqrt(gamma*d_out/(2*epsilon_0))*np.log(4*H/d_out)

def E_cj(Q):
    return np.sqrt(gamma/(epsilon_0*Rj(Q)))

def E_corner(V, R_corner, H):
    return V/(R_corner*np.log(4*H/R_corner))

def E_c(r):
    return 6.27e6*(1 + (8.8e-4/(r**(3/2))))

def Dd(Rj_val):
    return 1.89*2*Rj_val


def Score(E_cj_val, E_c_Rj_val, V_onset_val)



def filter_samples(samples):
    Q_min = np.sqrt(epsilon_r) * Q0

    d_in, d_out, R_corner = samples['d_in'], samples['d_out'], samples['R_corner']
    H, Q, V = samples['H'], samples['Q'], samples['V']

    mask = (
        (Q >= Q_min)
        & (V >= V_onset(d_out, H))
        & (E_cj(Q) < E_c(Rj(Q)))
        & (E_corner(V, R_corner, H) < E_c(R_corner))
        & (d_in < d_out)   # optional sanity check: inner diameter smaller than outer
    )
    return samples[mask].reset_index(drop=True)

# Scale the unit-cube samples to the physical bounds
param_names = list(boundary_parameters.keys())
l_bounds = np.array([boundary_parameters[k][0] for k in param_names])
u_bounds = np.array([boundary_parameters[k][1] for k in param_names])
steps    = np.array([resolution[k] for k in param_names])



# Number of distinguishable levels per parameter
n_levels = np.floor((u_bounds - l_bounds) / steps).astype(int) + 1
for k, n_k in zip(param_names, n_levels):
    print(f"{k:9s}: {n_k} levels")

# Latin Hypercube in the unit cube
sampler = qmc.LatinHypercube(d=len(param_names), optimization='random-cd', seed=42)
unit_samples = sampler.random(n=10000)

# Map each coordinate to an integer level index, then to a physical value on the grid
idx = np.minimum((unit_samples * n_levels).astype(int), n_levels - 1)
grid_samples = l_bounds + idx * steps

samples = pd.DataFrame(grid_samples, columns=param_names)

# Different unit-cube points can land on the same grid point, so remove duplicates
samples = samples.drop_duplicates().reset_index(drop=True)
print(f"{len(samples)} unique realizable designs")

filtered_samples = filter_samples(samples)
print(f"{len(filtered_samples)} satisfy the physical constraints")