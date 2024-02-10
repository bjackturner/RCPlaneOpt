#include <iostream>
#include<cmath>

/*
Air is a function that returns standard atmospheric conditions as a function of altitude (h).
This function uses an atmospheric model from NASA and is accurate up to 35 km, which is too high for most 
planes. This model will be used for all subsequent air properties in future calculations. The program will 
return air pressure (P), density (rho), temperature (T), dynamic viscosity (mu), the ratio of specific 
heat (gamma), gas constant (R), and acceleration of gravity (g). Like the rest of the program, all values 
must be calculated in standard metric units. 
*/

double Air(double h) {

    double T; double P; double rho; double mu;

    double R = 287.052874247; // Gas Constant (Air)
    double gamma = 1.4; // Specific Heat Ratio

    double g = pow(((6.6743e-11) * (5.9722e24))/((6.371000e6) + h),2); // Gravitational Acceleration

    if (h <= 11000) { // Troposphere
        T = 288.15 - 0.00649 * h; // Temperature Model
        P = pow(101325 * (T / 288.15),5.256); // Pressure Model
        rho = P/(R * T); // Density model
        mu = 0.00001789 * (pow((T / 288.15),1.5)) * (398.55 / (T + 110.4)); // Viscosity Model
    }

    else if (h > 11000 && h <= 25000) { // Lower Stratosphere
        T = 216.69; // Temperature Model
        P = 22632.06 * exp(1.73 - 0.000157 * h); // Pressure Model
        rho = P/(R * T); // Density model
        mu = 0.00001789 * (pow((T / 288.15),1.5)) * (398.55 / (T + 110.4)); // Viscosity Model
    }

    else if (h > 25000) { // Upper Stratosphere
        T = 141.94 + 0.00299 * h; // Temperature Model
        P = pow(2488 * (T / 216.6),-11.388); // Pressure Model
        rho = P/(R * T); // Density model
        mu = 0.00001789 * (pow((T / 288.15),1.5)) * (398.55 / (T + 110.4)); // Viscosity Model
    };

    return P, rho, T, mu, gamma, R, g;
}
