import numpy as np
import matplotlib.pyplot as plt

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
    position += (dt / 6.0) * (velocity + 2.0 * velocity_k2 + 2.0 * velocity_k3 + velocity_k4)
    velocity += (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return position, velocity

# 
def euler_integration(position, velocity, mass, airDensity, g, S, muRolling, alpha, aeroData, thrustData, dt):
    # Calculate force at the current position and velocity
    force = calculateForce(mass, airDensity, g, S, muRolling, alpha, position, velocity, aeroData, thrustData)
    
    # Update velocity using Euler's method
    new_velocity = velocity + (force / mass) * dt
    
    # Update position using the updated velocity
    new_position = position + new_velocity * dt
    
    return new_position, new_velocity

# Initialize kinematic state vectors 
def initializeVectors():
    
    t = 0 # Time vector (s)
    position = np.array([0.0,0.0]) # Position vector (X,Y)
    velocity = np.array([1e-10,0.0]) # Velocity vector (X,Y)
    acceleration = np.array([0.0,0.0]) # Acceleration vector (X,Y)

    return t, position, velocity, acceleration

# Ensure body does not travel through floor
def enforceFloorBoundary(position, velocity, floorHeight):

    # Check if position is below floor y < 0 
    if position[1] < floorHeight:
        position[1] = 0 # Reset position on floor
        velocity[1] = abs(velocity[1]) * 0.4 # Invert vertical velocity and apply collision damping

    return position, velocity

# Apply friction force as long as object is on ground and moving
def applyFriction(position, velocity, force, muRolling):
    
    # Apply friction force in the x-direction (object moving forward)
    if position[1] <= 0 and velocity[0] > 0:
        force[0] -= muRolling * abs(force[1]) # Apply friction force backward

    # Apply friction force in the x-direction (object moving backward)
    elif position[1] <= 0 and velocity[0] < 0:
        force[0] += muRolling * abs(force[1]) # Apply friction force forward

    return force

#
def calculateAero(aeroData, velocity, alpha):
    Cl = 1.2
    Cd = 0.05

    return Cl, Cd

#
def calculateThrust(thrustData, velocity, theta, alpha):

    effVelocity = np.dot(np.array([[np.cos(alpha+theta), -np.sin(alpha+theta)],[np.sin(alpha+theta), np.cos(alpha+theta)]]), np.transpose(velocity))

    return 20

#
def calculateForce(mass, airDensity, g, S, muRolling, alpha, position, velocity, aeroData, thrustData):

    theta = np.arctan(velocity[1]/velocity[0]) # Velocity vector angle (from horizontal)

    Cl, Cd = calculateAero(aeroData, velocity, alpha) #
    thrust = calculateThrust(thrustData, velocity, theta, alpha) # 

    force = np.array([0,0]) # Regenerate empty force array (X,Y)

    force[0] = (1/2 * airDensity * (np.linalg.norm(velocity))**2 * S) * (-Cl * np.sin(theta) - Cd * np.cos(theta)) + (thrust * np.cos(alpha + theta)) # Force in X direction (N)
    force[1] = (1/2 * airDensity * (np.linalg.norm(velocity))**2 * S) * (Cl * np.cos(theta) - Cd * np.sin(theta)) - (mass * g) + (thrust * np.sin(alpha + theta)) # Force in Y direction (N)

    force = applyFriction(position, velocity, force, muRolling) # Apply friction if aircraft is on the ground

    return force

# Takeoff distance numerical intergration
def takeoffDistance(geometricData, aeroData, thrustData):


    t, position, velocity, acceleration = initializeVectors()

    # Initialize additional state vectors
    force = np.array([0,0,0])

    print(position, velocity, acceleration)



    groundRoll = 1

    return groundRoll


t, position, velocity, acceleration = initializeVectors()

# Arrays to store position and time for plotting
positions = [position.copy()]
times = [0.0]
floorHeight = 0

# Simulation loop
t = 0.0
dt = 0.02

while t < 10:

    position, velocity = euler_integration(position, velocity, 10, 1.225, 9.81, 2, 0.05, 0, 2, 2, dt)
    # position, velocity = rk4_integration(position, velocity, 10, 1.225, 9.81, 2, 0.05, 0.05, 2, 2, dt)

    position, velocity = enforceFloorBoundary(position, velocity, floorHeight)

    # Append current position and time to arrays
    positions.append(position.copy())
    times.append(t)

    # Update time
    t += dt
    print(t)

# Convert lists to NumPy arrays for plotting
positions = np.array(positions)
times = np.array(times)

# Plotting
plt.plot(positions[:, 0], positions[:, 1])
plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.title('Trajectory of Object Falling and Bouncing Off Floor with Gravity')
plt.grid(True)
plt.show()