#include <iostream>
#include <cmath>
#include <algorithm>

double AirDensity(double h);
double Gravity(double h);

double mag(double x[]);

double Takeoff(double h, double m, double AeroData, double ThrustData) {

    double rho = AirDensity(h);
    double g = Gravity(h); // Gravitational Acceleration (m\s^2)

    double dt = 10e-3; // Time Step (s)
    double F[2] = {0,0}; // Set Initial Force Vector (N)

    double x[2] = {0,0}; // Set Initial Position Vector (m)
    double v[2] = {0,0}; // Set Initial Velocity Vector (m/s)
    double a[2] = {0,0}; // Set Initial Acceleration Vector (m/s^2)
    
    while(x[1] <= 0) {

        F[1] = (1/2 * rho * mag(v)^2 * CL * S) * cos(atan(v[1]/v[0])) - (m * g) - (1/2 * rho * mag(v)^2 * CD * S)*sin(atan(v[1]/v[0])) + (Thrust * sin(alpha + atan(v[1]/v[0]))); // Force Y direct (N)
        F[0] = -(1/2 * rho * mag(v)^2 * CL * S) * sin(atan(v[1]/v[0])) - (1/2 * rho * mag(v)^2 * CD * S)*cos(atan(v[1]/v[0])) + (Thrust * cos(alpha + atan(v[1]/v[0]))); // Force X direct (N)
       
        if (x[1] <= 0) {
            F[0] += std::min(F[1] * muRoll, 0);
            F[1] += std::max(F[1], 0);

        }

        for(int i = 0; i < 2; i++) {
            a[i] = F[i]/m; // Update Acceleration (m/s^2)
            v[i] = v[i] + a[i]*dt; // Update Velocity (m/s)
            x[i] = x[i] + v[i]*dt; // Update Position (m)

        }

    }

    return x[0]; // Return Ground Roll
}
