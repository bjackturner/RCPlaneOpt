import numpy as np

"""
Flow Properties contain functions and various flow properties in subsequent calculations. 
The functions included are an atmospheric model from NASA, Mach, and Reynolds number calculators. 
"""

# Static temperature (k) model as a function of altitude (m) using standard atm
def airTemp(h):
    # Troposphere properties
    if h <= 11000:  
        return 288.15 - 0.00649 * h 

    # Lower stratosphere properties
    elif 11000 < h <= 25000: 
        return 216.69 

    # Upper stratosphere properties
    else:  
        return 141.94 + 0.00299 * h 

# Static pressure (Pa) model as a function of altitude (m) using standard atm
def airPressure(h):

    T = airTemp(h) # Get air temp

    # Troposphere properties
    if h <= 11000: 
        return 101325 * (T / 288.15)**5.256

    # Lower stratosphere properties
    elif 11000 < h <= 25000: 
        return 22632.06 * np.exp(1.73 - 0.000157 * h)

    # Upper stratosphere properties
    else:  
        return 2488 * (T / 216.6)**-11.388
    
# Static density (kg/m^3) model as a function of altitude (m) using standard atm
def airDensity(h):

    T = airTemp(h) # Get air temp
    P = airPressure(h) # Get air pressure
    R = 287.052874247 # Gas constant (air)

    return P / (R * T)

# Dynamic viscosity (Pa/s) model as a function of altitude (m) using standard atm
def airDynamicViscosity(h):

    T = airTemp(h) # Get air temp

    return 0.00001789 * (T / 288.15)**1.5 * (398.55 / (T + 110.4))

# Gravitational acceleration (m/s^2) model as a function of altitude (m) using standard atm
def gravity(h):

    return (6.6743e-11 * 5.9722e24) / (6.371000e6 + h)**2

# Calculates Reynolds number ((rho * u * l) / mu)
def reynoldsNum(u, l, rho, mu):

    return (rho * u * l) / mu  # Reynolds number calculation

# Calculates Mach number (u / a)
def machNum(u, gamma, R, T):

    a = np.sqrt(gamma * R * T)  # Speed of Sound calculation

    return u / a # Mach Number calculation