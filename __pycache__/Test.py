import numpy as np

# Given points and values
x1 = 0
y1 = 0
m1 = 0
x2 = 1
y2 = 1
m2 = 2
y1_prime = -(x1) + 1
y2_prime = -(x2) + 1

# Calculate the slope m of the linear line
m = (y2_prime - y1_prime) / (x2 - x1)


# Calculate the y-intercept c of the linear line
c = y1_prime - m * x1

# Solve for x-coordinate of the intersection point by setting the two functions equal to each other
x_intersection = (c - y1 + m1 * x1) / (m1 - m)

# Calculate y-coordinate of the intersection point using either of the functions
y_intersection = m * x_intersection + c

# Print the intersection point
print("Intersection point:", x_intersection, y_intersection)