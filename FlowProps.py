import numpy as np

"""
Flow Properties contain functions and various flow properties in subsequent calculations. 
The functions included are an atmospheric model from NASA, Mach, and Reynolds number calculators. 
"""

# Returns properties of air
def air(h):
    
    R = 287.052874247 # Gas Constant (Air)
    gamma = 1.4 # Specific Heat Ratio

    g = (6.6743e-11 * 5.9722e24) / (6.371000e6 + h)**2 # Gravitational Acceleration

    # Troposphere Properties
    if h <= 11000:  
        T = 288.15 - 0.00649 * h # Temperature Model
        P = 101325 * (T / 288.15)**5.256 # Pressure Model
        rho = P / (R * T) # Density model
        mu = 0.00001789 * (T / 288.15)**1.5 * (398.55 / (T + 110.4)) # Viscosity Model

    # Lower Stratosphere Properties
    elif 11000 < h <= 25000: 
        T = 216.69 # Temperature Model
        P = 22632.06 * np.exp(1.73 - 0.000157 * h) # Pressure Model
        rho = P / (R * T) # Density model
        mu = 0.00001789 * (T / 288.15)**1.5 * (398.55 / (T + 110.4)) # Viscosity Model

    # Upper Stratosphere Properties
    else:  
        T = 141.94 + 0.00299 * h # Temperature Model
        P = 2488 * (T / 216.6)**-11.388 # Pressure Model
        rho = P / (R * T) # Density model
        mu = 0.00001789 * (T / 288.15)**1.5 * (398.55 / (T + 110.4)) # Viscosity Model

    return P, rho, T, mu, gamma, R, g # Return Air Properties


# Calculates Reynolds Number.
def reynolds(u, l, rho, mu):

    return (rho * u * l) / mu  # Reynolds number calculation


# Calculates Mach Number.
def mach(u, gamma, R, T):

    a = np.sqrt(gamma * R * T)  # Speed of Sound calculation

    return u / a # Mach Number calculation