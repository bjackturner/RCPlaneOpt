import numpy as np

alpha = 10 /180 * np.pi

point1 = [1,0]
point2 = (0.5, -1e-7)

dx = point2[0] - point1[0]
dy = point2[1] - point1[1]

theta = np.arctan2(dy,dx) + np.pi/2 - alpha

if theta < 0:
    theta += 2*np.pi

print(theta*180/np.pi)