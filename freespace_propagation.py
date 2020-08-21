# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 16:52:05 2020

@author: Benjamin
"""
"""
    Propagation through free space. Refractive index taken to be one.
"""


import numpy as np
import matplotlib.pyplot as plt

def waist_angle_from_q_f(q_f, index, wavelength):
    # Calculate new waist and angle from q_f
    q_f_inverse = 1 / q_f
    radius = 1 / q_f_inverse.real 
    waist = np.sqrt(- wavelength / (q_f_inverse.imag * np.pi * index))
    angle = np.arctan(waist / radius)
    return waist, angle
    
    
def q_i_from_waist(waist, angle, index, wavelength):
    # Calculate q_f from waist
    radius = waist / np.tan(angle) # Radius of curvature using sohcahtoa
    q_i_inverse = 1 / radius - 1j * (wavelength)/(index * np.pi * waist**2)
    return 1 / q_i_inverse


def q_f_from_matrix(q_i, matrix):
    # Extracts matrix elements and calculates q_f
    A = matrix[0, 0]
    B = matrix[0, 1]
    C = matrix[1, 0]
    D = matrix[1, 1]
    q_f = (A * q_i + B) / (C * q_i + D)
    return q_f

# Define initial variables 
wavelength = 8e-7
waist_0 = 180e-6 # initial beam waist
angle_0 = 0 # initial half apex angle of beam
index = 1 # refractive index of free space
theta = 0 # half apex angle of beam

# Define distance range over 300mm total
z = np.linspace(start = 0, stop = 0.3, num = 1001)
# Ensure all values are to 4dp
z = np.round(z, 4)

# Calculate rayleigh length
z_r = (np.pi * index * (waist_0)**2) / wavelength

# Initial complex beam parameter
q_i = z[0] + z_r * 1j

# Take distance to be the last value within z
distance = z[-1]

# Define beam waist to be an empty array with same length as z
waist = np.zeros(len(z))
# The first beam waist must be waist_0
waist[0] = waist_0

# Define beam waist to be an empty array with same length as z
angle = np.zeros(len(z))
# The first beam waist must be waist_0
angle[0] = angle_0

# Define step to be z[1] - z[0] (arbitrary, all steps are uniform)
step = z[1] - z[0]

# Define matrix of system according to ray transfer matrix analysis
matrix = np.array([[1, step], [0, 1]]) # ABCD matrix for free space prop.

matrix_init = np.identity(2)

# Iterate through waist to calculate its value in each z
for i in range(len(z)):
        
    # Dont need to calculate new parameters if at end of loop
    if i == max(range(len(z))):
        continue
        

    # update current position, waist, angle
    current_z = z[i]
    current_waist = waist[i]
    current_angle = angle[i]
    
    # calculate new matrix
    if i == 0:
        current_matrix = np.dot(matrix, matrix_init)
    else:
        current_matrix = np.dot(matrix, current_matrix)

    
    # q_i is based on current position
    q_i = q_i_from_waist(waist = waist_0, 
                       angle = angle_0,
                       index = index,
                       wavelength = wavelength)

    q_f = q_f_from_matrix(q_i, current_matrix)

    new_waist, new_angle = waist_angle_from_q_f(q_f, index, wavelength)

    # Write new waist and angle
    waist[i+1] = new_waist
    angle[i+1] = new_angle

    
plt.plot(z, np.zeros(len(z)), 'k--')
plt.plot(z, waist, 'b-')
plt.plot(z, -waist, 'b-')
plt.fill_between(z, waist, -waist, facecolor='blue', alpha=0.2)

plt.show()
    