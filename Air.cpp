#include <iostream>
#include<cmath>

/*
Air includes functions that returns standard atmospheric conditions as a function of altitude (h).
The functions use an atmospheric model from NASA, which accurate up to 35 km. This model will be used 
for all subsequent air properties in future calculations. The program will return air pressure (P), 
density (rho), temperature (T), dynamic viscosity (mu), the ratio of specific heats (gamma), gas 
constant (R), and acceleration of gravity (g). Like the rest of the program, all values 
must be calculated in standard metric units. 
*/

// Air Temperature as a function of altitude  
double AirTemp(double h) {

    if (h <= 11000) { // Troposphere
        return 288.15 - 0.00649 * h; // Temperature Model
    }

    else if (h > 11000 && h <= 25000) { // Lower Stratosphere
        return 216.69; // Temperature Model
    }

    else if (h > 25000) { // Upper Stratosphere
        return 141.94 + 0.00299 * h; // Temperature Model
    }
}

// Air Pressure as a function of altitude  
double AirPressure(double h) {

    double T = AirTemp(h); // Find Air Temperature (func of Alt)

    if (h <= 11000) { // Troposphere
        return 101325 * pow((T / 288.15), 5.256); // Pressure Model
    }

    else if (h > 11000 && h <= 25000) { // Lower Stratosphere
        return 22632.06 * exp(1.73 - 0.000157 * h); // Pressure Model
    }

    else if (h > 25000) { // Upper Stratosphere
        return 2488 * pow((T / 216.6), -11.388); // Pressure Model
    }
}

// Air Density as a function of altitude  
double AirDensity(double h) {

    double T = AirTemp(h); // Find Air Temperature (func of Alt)
    double P = AirPressure(h); // Find Air Pressure (func of Alt)

    return P/(287.052874247 * T); // Density model (Ideal Gas)
}

// Air Viscosity as a function of altitude  
double AirViscosity(double h) {

    double T = AirTemp(h); // Find Air Temperature (func of Alt)

    if (h <= 11000) { // Troposphere
        return 0.00001789 * (pow((T / 288.15),1.5)) * (398.55 / (T + 110.4)); // Viscosity Model
    }

    else if (h > 11000 && h <= 25000) { // Lower Stratosphere
        return 0.00001789 * (pow((T / 288.15),1.5)) * (398.55 / (T + 110.4)); // Viscosity Model
    }

    else if (h > 25000) { // Upper Stratosphere
        return 0.00001789 * (pow((T / 288.15),1.5)) * (398.55 / (T + 110.4)); // Viscosity Model
    }
}

// Gravitational acceleration as a function of altitude  
double Gravity(double h) {

    return ((6.6743e-11) * (5.9722e24))/pow(((6.371000e6) + h),2); // Gravitational Acceleration
}
