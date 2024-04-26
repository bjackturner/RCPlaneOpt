import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as sci
import time

import FlowProps
import PropCalc
import Preprocessing

# Create 4, 5, and 6 series NACA airfoils (NACA digits [string], total number of points)
def buildNACAFoil(numNACA,numPoints,closedTip):

    NACAParmaters = [] # Initialize NACA Parameters array

    chord = np.linspace(-np.pi,np.pi,numPoints) # Create array of points from -pi to pi for cosine distribution
    chord = np.sign(chord) * 0.5*(1 - np.cos(chord)) # Finilaize cosine distribution from -1 to 1
    
    airfoil = np.zeros((2, numPoints)) # Initialize airfoil array

    for char in numNACA:
        NACAParmaters.append(int(char)) # Determine how many inputs given
    
    # Four digit NACA airfoils
    if len(NACAParmaters) == 4:

        # Thickness distribution constants
        a0 = 0.2969
        a1 = -0.1260
        a2 = -0.3516
        a3 = 0.2843
        
        if closedTip == True:
            a4 = -0.1036

        elif closedTip == False:
            a4 = -0.1015

        i = 0 # Initialize counting variable :(

        # Extract airfoil paramaters from input
        maxCamber = int(numNACA[0])/100 # Max camber % of chord
        percentMaxCamber = int(NACAParmaters[1])/10 # Location of max camber % of chord
        thickness = (10*int(NACAParmaters[2]) + int(NACAParmaters[3]))/100 # Thickness % of chord

        # Generate four digit airfoil (loop through each x coordinate)
        for x in chord:
            # If x coordinate is ahead of location of maximum camber (ensure first two NACA digits are not zero)
            if np.abs(x) <= percentMaxCamber and maxCamber*percentMaxCamber != 0:
                yc = maxCamber/percentMaxCamber**2 * (2*percentMaxCamber*np.abs(x) - np.abs(x)**2) # Camber line as function of x
                theta = np.arctan(((2*maxCamber)/percentMaxCamber**2) * (percentMaxCamber - np.abs(x))) # Angle of camber line (derivative of camber line function)

            # If x coordinate is behind location of maximum camber (ensure first two NACA digits are not zero)
            elif np.abs(x) > percentMaxCamber and maxCamber*percentMaxCamber != 0:
                yc = maxCamber/(1 - percentMaxCamber)**2 * ((1 - 2*percentMaxCamber) + 2*percentMaxCamber*np.abs(x) - np.abs(x)**2) # Camber line as function of x
                theta = np.arctan(((2*maxCamber)/(1-percentMaxCamber)**2) * (percentMaxCamber - np.abs(x))) # Angle of camber line (derivative of camber line function)

            # If airfoil is symmetric (ensure first two NACA digits are zero)
            elif maxCamber == 0 and percentMaxCamber ==0:
                yc = 0 # Camber line as function of x 
                theta = 0 # Derivative of camber line

            # Error message (stops code in the case of invalid inputs)
            else:
                print('Error. Invalid Input') # Print error message
                break # Break loop and end program 

            # Generate top side of airfoil (trailing edge to leading edge)
            if x >= 0:
                x = np.abs(x) # Ensure x-value is positive
                yt = 5*thickness * (a0*np.sqrt(x) + a1*x + a2*x**2 + a3*x**3 + a4*x**4) # Thickness as a function x
                
                airfoil[0][i] = x - yt*np.sin(theta) # Airfoil x coordinate 
                airfoil[1][i] = yc + yt*np.cos(theta) # Airfoil y coordinate

            # Generate bottom side of airfoil (leading edge to trailing edge)
            elif x < 0:
                x = np.abs(x) # Ensure x-value is positive
                yt = 5*thickness * (a0*np.sqrt(x) + a1*x + a2*x**2 + a3*x**3 + a4*x**4) # Thickness as a function x

                airfoil[0][i] = x + yt*np.sin(theta) # Airfoil x coordinate
                airfoil[1][i] = yc - yt*np.cos(theta) # Airfoil y coordinate

            i += 1 # Counting variable :(
            
        return airfoil # Return 4 digit NACA

    # Five digit NACA airfoils
    elif len(NACAParmaters) == 5:

        # Thickness distribution constants
        a0 = 0.2969
        a1 = -0.1260
        a2 = -0.3516
        a3 = 0.2843
        
        if closedTip == True:
            a4 = -0.1036

        elif closedTip == False:
            a4 = -0.1015

        i = 0 # Initialize counting variable :(

        # Extract airfoil paramaters from input
        designCL = int(numNACA[0]) * 3 / 20 # Design lift coefficient
        percentMaxCamber = int(NACAParmaters[1])/20 # Max camber % of chord
        reflex = int(NACAParmaters[2]) # Reflex Camber [Reflex = 0 - Standard camber | Reflex = 1 - reflexed camber]
        thickness = (10*int(NACAParmaters[3]) + int(NACAParmaters[4]))/100 # Thickness % of chord

        #   
        if reflex == 0:

            def solveR(m):
                return (m*(1 - np.sqrt(m/3))) - percentMaxCamber

            m = sci.fsolve(solveR, 0)
            N = ((3*m - 7*m**2 + 8*m**3 - 4*m**4)/(np.sqrt(m - m**2))) - (3/2 * (1 - 2*m) * (np.pi/2 - np.arcsin(1 - 2*m)))
            k1 = 6 * designCL / N

            for x in chord:  
                # If x coordinate is ahead of location of maximum camber (ensure first two NACA digits are not zero)
                if np.abs(x) < m and percentMaxCamber*designCL != 0:
                    yc = k1/6 * (np.abs(x)**3 - 3*m*np.abs(x)**2 + m**2*(3 - m)*np.abs(x)) # Camber line as function of x
                    theta = np.arctan(k1/6 * (3*np.abs(x)**2 - 6*m*np.abs(x) + (3-m)*m**2)) # Angle of camber line (derivative of camber line function)

                # If x coordinate is behind location of maximum camber (ensure first two NACA digits are not zero)
                elif np.abs(x) >= m and percentMaxCamber*designCL != 0:
                    yc = k1*m**3 / 6 * (1 - np.abs(x)) # Camber line as function of x
                    theta = np.arctan(-k1 * m**3 / 6) # Angle of camber line (derivative of camber line function)

                # Error message (stops code in the case of invalid inputs)
                else:
                    print('Error. Invalid Input') # Print error message
                    break # Break loop and end program 

                # Generate top side of airfoil (trailing edge to leading edge)
                if x >= 0:
                    x = np.abs(x) # Ensure x-value is positive
                    yt = 5*thickness * (a0*np.sqrt(x) + a1*x + a2*x**2 + a3*x**3 + a4*x**4) # Thickness as a function x
                    
                    airfoil[0][i] = x - yt*np.sin(theta) # Airfoil x coordinate 
                    airfoil[1][i] = yc + yt*np.cos(theta) # Airfoil y coordinate

                # Generate bottom side of airfoil (leading edge to trailing edge)
                elif x < 0:
                    x = np.abs(x) # Ensure x-value is positive
                    yt = 5*thickness * (a0*np.sqrt(x) + a1*x + a2*x**2 + a3*x**3 + a4*x**4) # Thickness as a function x

                    airfoil[0][i] = x + yt*np.sin(theta) # Airfoil x coordinate
                    airfoil[1][i] = yc - yt*np.cos(theta) # Airfoil y coordinate

                i += 1 # Counting variable :(

        elif reflex == 1:

            def solveR(m):
                return (m*(1 - np.sqrt(m/3))) - percentMaxCamber

            m = sci.fsolve(solveR, 0)
            N = ((3*m - 7*m**2 + 8*m**3 - 4*m**4)/(np.sqrt(m - m**2))) - (3/2 * (1 - 2*m) * (np.pi/2 - np.arcsin(1 - 2*m)))
            k1 = 6 * designCL / N
            k21 = (3 * (m - percentMaxCamber)**2 - m**3) / (1 - m)**3

            for x in chord:
                if np.abs(x) < m and percentMaxCamber*designCL != 0:
                    yc = k1/6 * ((np.abs(x) - m)**3 - np.abs(x)*k21*(1 - m)**3 - np.abs(x)*m**3 + m**3) # Camber line as function of x
                    theta = np.arctan(k1/6 * (3*(np.abs(x) - m)**2 - k21*(1 - m)**3 - m**3)) # Angle of camber line (derivative of camber line function)

                # If x coordinate is behind location of maximum camber (ensure first two NACA digits are not zero)
                elif np.abs(x) >= m and percentMaxCamber*designCL != 0:
                    yc = k1/6 * (k21*(np.abs(x) - m)**3 - np.abs(x)*k21*(1 - m)**3 - np.abs(x)*m**3 + m**3) # Camber line as function of x
                    theta = np.arctan(k1/6 * (3*k21*(np.abs(x) - m)**2 - k21*(1 - m)**3 - m**3)) # Angle of camber line (derivative of camber line function)

                # Error message (stops code in the case of invalid inputs)
                else:
                    print('Error. Invalid Input') # Print error message
                    break # Break loop and end program 

                # Generate top side of airfoil (trailing edge to leading edge)
                if x >= 0:
                    x = np.abs(x) # Ensure x-value is positive
                    yt = 5*thickness * (a0*np.sqrt(x) + a1*x + a2*x**2 + a3*x**3 + a4*x**4) # Thickness as a function x
                    
                    airfoil[0][i] = x - yt*np.sin(theta) # Airfoil x coordinate 
                    airfoil[1][i] = yc + yt*np.cos(theta) # Airfoil y coordinate

                # Generate bottom side of airfoil (leading edge to trailing edge)
                elif x < 0:
                    x = np.abs(x) # Ensure x-value is positive
                    yt = 5*thickness * (a0*np.sqrt(x) + a1*x + a2*x**2 + a3*x**3 + a4*x**4) # Thickness as a function x

                    airfoil[0][i] = x + yt*np.sin(theta) # Airfoil x coordinate
                    airfoil[1][i] = yc - yt*np.cos(theta) # Airfoil y coordinate

                i += 1 # Counting variable :(
                
        else:
            print('Error. Invalid Input') # Print error message
            return # Break loop and end program 

        return airfoil

    # Six digit NACA airfoils
    elif len(NACAParmaters) == 6:
        print('6 digit Foil')

    # Invalid input statement
    else:
        print('Error. Invalid Input')

def buildPannelGeomerty(airfoil,alpha):

    numPannles = len(airfoil[1]) - 1

    S = np.zeros(numPannles)
    X = np.zeros(numPannles)
    Y = np.zeros(numPannles)
    phi = np.zeros(numPannles)

    for i in range(numPannles):
        X[i] = (airfoil[0,i+1] + airfoil[0,i])/2
        Y[i] = (airfoil[1,i+1] + airfoil[1,i])/2

        dx = airfoil[0,i+1] - airfoil[0,i]
        dy = airfoil[1,i+1] - airfoil[1,i]

        S[i] = np.sqrt((dx)**2 + (dy)**2)

        phi[i] = np.arctan2(dy,dx) + np.pi/2 - alpha

        if phi[i] < 0:
            phi[i] += 2*np.pi

    return X, Y, S, phi


# Invert airfoil (airfoil)
def invertAirfoil(airfoil):

    airfoil[1,:] = airfoil[1,:] * -1 # Inverted airfoil by taking the opposite y coordinate

    return airfoil # Return inverted airfoil

# Rotate airfoil (airfoil, rotation angle [rads], rotation location [% of chord])
def rotateAirfoil(airfoil,rotationAngle,rotationLocation=None):

    # currentRotation = -np.arctan((airfoil[1, int(len((airfoil[1]))/2)] - airfoil[1, 0]) / (airfoil[0, int(len((airfoil[1]) + 1)/2)] - airfoil[0, 0])) * 0# Find chord length to correct airfoil scale for rotation location
    # print('Current Rotation Angle = ', currentRotation * 180/ np.pi)

    # chordLength = np.sqrt((airfoil[0, int(len(airfoil[1])/2)] - airfoil[0, 0])**2 + (airfoil[1, int(len((airfoil[1]))/2)] - airfoil[1, 0])**2) # Find chord length to correct airfoil scale for rotation location
    # print('chord for Rotation = ',chordLength)
    # print('-')

    chordLength = 1
    currentRotation = 0

    # Build rotation matrix
    R = np.array([[np.cos(-rotationAngle + currentRotation), -np.sin(-rotationAngle + currentRotation)],
                  [np.sin(-rotationAngle + currentRotation), np.cos(-rotationAngle + currentRotation)]])
    
    airfoil[0, :] -= rotationLocation * chordLength
    # airfoil[1, :] -= rotationLocation * chordLength

    airfoil = np.dot(R,airfoil)

    airfoil[0, :] += rotationLocation * chordLength 
    # airfoil[1, :] += rotationLocation * chordLength

    return airfoil

# Generate flap on pre-existing airfoil (airfoil, flap location [% of chord], flap angle [degrees])
def flapAirfoil(airfoil,flapLocation,flapAngle):

    return airfoil

# Scale airfoil (airfoil, scale factor)
def scaleAirfoil(airfoil,scaleFactor):

    airfoil[0,:]  = airfoil[0,:] * scaleFactor
    airfoil[1,:]  = airfoil[1,:] * scaleFactor

    return airfoil # Return scaled airfoil

airfoil = buildNACAFoil('24012',31,True)
# test = invertAirfoil(airfoil)
# test = rotateAirfoil(airfoil, 10/180*np.pi, 0.25)
# test = scaleAirfoil(test, 2)

X, Y, S, phi = buildPannelGeomerty(airfoil,0)


# with open("RCPlaneOpt/n23012.txt", "r") as file:
#     lines = file.readlines()

# # Extract x and y coordinates from the data
# x = []
# y = []
# for line in lines:
#     parts = line.split()
#     x.append(float(parts[0]))
#     y.append(float(parts[1]))

# plt.plot(x,y)
plt.plot(airfoil[0,:],airfoil[1,:])

XC = np.zeros(2)
YC = np.zeros(2)
for i in range(len(airfoil[1]) - 1):
    XC[0] = X[i]
    XC[1] = X[i] + S[i]*np.cos(phi[i])
    YC[0] = Y[i]
    YC[1] = Y[i] + S[i]*np.sin(phi[i])
    if (i == 0):
        plt.plot(XC,YC,'b-',label='First Panel')
    elif (i == 1):
        plt.plot(XC,YC,'g-',label='Second Panel')
    else:
        plt.plot(XC,YC,'r-')

plt.xlabel('X-Axis')
plt.ylabel('Y-Axis')
plt.title('Panel Geometry')
plt.axis('equal')
plt.show()                   
