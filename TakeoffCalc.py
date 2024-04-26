import numpy as np
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation

import FlowProps
import PropCalc
import Preprocessing

# 
def rk4_integration(position, velocity, mass, airDensity, g, S, muRolling, alpha, aeroData, thrustData, dt):
    # Calculate forces at the initial position and velocity
    k1 = calculateForce(mass, airDensity, g, S, muRolling, alpha, position, velocity, aeroData, thrustData)
    
    # Use k1 to estimate k2 at the midpoint of the time interval
    position_k2 = position + 0.5 * dt * velocity
    velocity_k2 = velocity + 0.5 * dt * k1
    k2 = calculateForce(mass, airDensity, g, S, muRolling, alpha, position_k2, velocity_k2, aeroData, thrustData)
    
    # Use k2 to estimate k3 at the midpoint of the time interval
    position_k3 = position + 0.5 * dt * velocity_k2
    velocity_k3 = velocity + 0.5 * dt * k2
    k3 = calculateForce(mass, airDensity, g, S, muRolling, alpha, position_k3, velocity_k3, aeroData, thrustData)
    
    # Use k3 to estimate k4 at the end of the time interval
    position_k4 = position + dt * velocity_k3
    velocity_k4 = velocity + dt * k3
    k4 = calculateForce(mass, airDensity, g, S, muRolling, alpha, position_k4, velocity_k4, aeroData, thrustData)

    # Update position and velocity using the RK4 method
    positionDt = (dt / 6.0) * (velocity + 2.0 * velocity_k2 + 2.0 * velocity_k3 + velocity_k4)
    velocityDt = (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    position, velocity = enforceFloorBoundary(position, velocity, positionDt, velocityDt, dt)

    return position, velocity

# 
def euler_integration(position, velocity, mass, airDensity, g, S, muRolling, alpha, aeroData, thrustData, dt):
    # Calculate force at the current position and velocity
    force = calculateForce(mass, airDensity, g, S, muRolling, alpha, position, velocity, aeroData, thrustData, dt)
    
    # Update velocity using Euler's method
    velocityDt = velocity + (force / mass) * dt
    
    # Update position using the updated velocity
    positionDt = position + velocityDt * dt

    position, velocity = enforceFloorBoundary(position, positionDt, velocityDt, dt)


    return position, velocity

# Used to find the unit vector normal to a line
def calculateBounceAngle(point1, point2):

    theta = np.arctan((point2[1] - point1[1])/(point2[0] - point1[0]))

    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])

# Initialize kinematic state vectors 
def initializeVectors():
    
    t = 0 # Time vector (s)
    position = np.array([0.0,0.0]) # Position vector (X,Y)
    velocity = np.array([0.0,0.0]) # Velocity vector (X,Y)
    acceleration = np.array([0.0,0.0]) # Acceleration vector (X,Y)

    return t, position, velocity, acceleration

def calculateFloorHeight(position):

    return np.sin( position/10)

# Ensure body does not travel through floorF
def enforceFloorBoundary(position, positionNext, velocityNext, dt):

    floorHeight = calculateFloorHeight(position[0]) # Find floor height at current position
    floorHeightNext = calculateFloorHeight(positionNext[0]) # Find floor height at next position

    if position[1] == floorHeight and positionNext[1] < floorHeightNext:
        # Reset position on floor
        positionNext[1] = floorHeightNext # Y floor position

        R = calculateBounceAngle([position[0], floorHeight], [positionNext[0], floorHeightNext])

        velocityNext = np.dot(velocityNext, R)

        velocityNext[1] = 0

        velocityNext = np.dot(R, velocityNext)

    # Check if next position is below floor at next time step
    elif positionNext[1] <= floorHeightNext:

        R = calculateBounceAngle([position[0], floorHeight], [positionNext[0], floorHeightNext])

        a1 = ((position[0] * positionNext[1]) - (position[1] * positionNext[0])) * (position[0] - positionNext[0])
        a2 = ((position[0] * floorHeightNext) - (floorHeight * positionNext[0])) * (position[0] - positionNext[0])
        a3 = ((position[0] - positionNext[0]) * (floorHeight - floorHeightNext)) - ((position[1] - positionNext[1]) * (position[0] - positionNext[0]))

        w1 = np.abs((position[0] - ((a1 - a2)/a3)) / (positionNext[0] - position[0]))
        w2 = np.abs((positionNext[0] - ((a1 - a2)/a3)) / (positionNext[0] - position[0]))

        # velocityNext[0] = ((w2 * velocity[0]) + (w1 * velocityNext[0]))/(w1 + w2)
        # velocityNext[1] = ((w2 * velocity[1]) + (w1 * velocityNext[1]))/(w1 + w2)

        positionNext[0] = ((a1 - a2)/a3)
        positionNext[1] = calculateFloorHeight(((a1 - a2)/a3))

        dt = ((w1 * 0) + (w2 * dt))/(w1 + w2)


        # print(velocityNext,'int')
        velocityNext = np.dot(velocityNext, R)

        velocityNext[1] = -velocityNext[1] * 0.8

        if velocityNext[1] < 1e-1:
            print(velocityNext)
            velocityNext = np.dot(R, velocityNext)
            

            positionNext = velocityNext * dt
            positionNext[1] = floorHeightNext

        else:
            print(velocityNext,'-')
            velocityNext = np.dot(R, velocityNext)
            positionNext += velocityNext * dt
            
    if positionNext[1] > floorHeightNext:
        print('takeoff = ',positionNext)
        
    return positionNext, velocityNext

def calculateDirectionVector(velocity):
    
    direction = np.array([0,0]) 

    direction[0] = np.dot(velocity,[1,0])
    direction[1] = np.dot(velocity,[0,1])

# Apply friction force as long as object is on ground and moving
def calculateFriction(position, velocity, force, muRolling, dt):

    if position[1] <= calculateFloorHeight(position[0]):
        
        floorHeight = calculateFloorHeight(position[0]) # Find floor height at current position
        R = calculateBounceAngle([position[0], floorHeight], [velocity[0] * dt - position[0], calculateFloorHeight(velocity[0] * dt - position[0])])

        force = np.dot(force, R)

        force[0] -= muRolling * abs(force[1]) * velocity[0]/np.abs(velocity[0])

        force = np.dot(R, force)

    return force

#
def calculateAero(aeroData, velocity, alpha):
    Cl = 1.2
    Cd = 0.05

    return Cl, Cd

#
def calculateThrust(thrustData, velocity, theta, alpha):

    effVelocity = np.dot(np.array([[np.cos(alpha+theta), -np.sin(alpha+theta)],[np.sin(alpha+theta), np.cos(alpha+theta)]]), np.transpose(velocity))

    return 0

#
def calculateForce(mass, airDensity, g, S, muRolling, alpha, position, velocity, aeroData, thrustData, dt):

    theta = np.arctan(velocity[1]/velocity[0])

    Cl, Cd = calculateAero(aeroData, velocity, alpha) #
    thrust = calculateThrust(thrustData, velocity, theta, alpha) # 

    cd = 0.05
    cl = 0.6

    force = np.array([0,0]) # Regenerate empty force array (X,Y)

    # force[0] = (1/2 * airDensity * (np.linalg.norm(velocity))**2 * S) * (-Cl * np.sin(theta) - Cd * np.cos(theta)) + (thrust * np.cos(alpha + theta)) # Force in X direction (N)
    # force[1] = (1/2 * airDensity * (np.linalg.norm(velocity))**2 * S) * (Cl * np.cos(theta) - Cd * np.sin(theta)) - (mass * g) + (thrust * np.sin(alpha + theta)) # Force in Y direction (N)

    force[0] = thrust + cd * airDensity * velocity[0]**2 / 2
    force[1] = -mass * g + cl * airDensity * velocity[0]**2 / 2

    force = calculateFriction(position, velocity, force, muRolling, dt) # Apply friction if aircraft is on the ground

    return force

# Takeoff distance numerical intergration
def takeoffDistance(geometricData, aeroData, thrustData):


    t, position, velocity, acceleration = initializeVectors()

    # Initialize additional state vectors
    force = np.array([0,0,0])


    groundRoll = 1

    return groundRoll


start_time = time.time()

t, position, velocity, acceleration = initializeVectors()

velocity[0] = 5
velocity[1] = 5
position[1] = 0

# Arrays to store position and time for plotting
positions = [position.copy()]
times = [0.0]

dt = 0.1

# while t < 30:

#     # print(position, velocity, t)

#     position, velocity = euler_integration(position, velocity, 10, 1.225, 9.81, 2, 0.05, 0, 2, 2, dt)
#     # position, velocity = rk4_integration(position, velocity, 10, 1.225, 9.81, 2, 0.05, 0.05, 2, 2, dt)

#     # Append current position and time to arrays
#     positions.append(position.copy())
#     times.append(t)

#     # Update time
#     t += dt



# Function to update animation frames
def update(frame):
    global position, velocity, t

    # Update position and velocity using Euler integration
    position, velocity = euler_integration(position, velocity, 10, 1.225, 9.81, 2, 0.05, 0, 2, 2, dt)

    # Append current position and time to arrays
    positions.append(position.copy())
    times.append(t)

    
    # Update time
    t += dt

    # Update plot data
    line.set_data(position[0], position[1])

    return line,

def animate_trajectory():
    # Close existing figure to prevent the warning
    plt.close()

    # Create a figure and axis
    fig, ax = plt.subplots()
    ax.set_xlim(0, 200)
    ax.set_ylim(-3, 20)
    ax.set_xlabel('X Position (m)')
    ax.set_ylabel('Y Position (m)')
    ax.set_title('Trajectory of Object Falling and Bouncing Off Floor with Gravity')
    ax.grid(True)

    # Plot the floor
    x = np.linspace(-50, 200, 100)
    y = calculateFloorHeight(x)
    floor, = ax.plot(x, y)

    # Plot initial position
    global line
    line, = ax.plot(positions[0][0], positions[0][1], marker='o')

    # Create animation
    ani = FuncAnimation(fig, update, frames=10, blit=True)

    # Display animation
    plt.show()

# Call the function to display the animation
animate_trajectory()


# # Convert lists to NumPy arrays for plotting
# positions = np.array(positions)
# times = np.array(times)

# x = np.linspace(min(positions[:,0]),max(positions[:,0]),len(positions))
# y = calculateFloorHeight(x[:])

# end_time = time.time()
# execution_time = end_time - start_time

# print("Script execution time: {:.6f} seconds".format(execution_time))

# # Plotting
# plt.plot(positions[:, 0], positions[:, 1])
# plt.plot(x,y)
# plt.xlabel('X Position (m)')
# plt.ylabel('Y Position (m)')
# plt.title('Trajectory of Object Falling and Bouncing Off Floor with Gravity')
# plt.grid(True)
# plt.show()
