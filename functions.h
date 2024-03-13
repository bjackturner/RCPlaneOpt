#ifndef FUNCTIONS_H   
#define FUNCTIONS_H

#include "FlowProps.cpp"
#include "Air.cpp"

// Functions from FlowProps.cpp
double Mach(double u, double gamma, double R, double T);
double Reynolds(double rho, double u, double l, double mu);

// Functions from Air.cpp
double AirTemp(double h);
double AirPressure(double h);
double AirDensity(double h);
double AirViscosity(double h);
double Gravity(double h);

#endif