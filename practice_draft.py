# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 10:03:08 2022

@author: paloc
"""
import numpy as np
from math import pi
import matplotlib.pyplot as plt 
from Astronomy_functions import calculate_velocity_with_noise
Earth_mass = 6 * (10**24)

def calculate_amplitude(mass1, q_array , inclination_angle, period):
    """ Here is another function to calculate the amplitude of a binary system
    in a binary system.
    The inputs consist of:
        mass of star1 (bigger object)
        'q' which is the mass2:mass1 ratio np.random.uniform(0, 1, 6)
        angle of inclination
        and period of orbit in days, and the code will change it to seconds
        the amplitude output has units of km/sec
    """
    mass2 = q_array * mass1
    Gravitational_constant = 6.67431 * (10**-11)
    orbital_period = period * 86400  #change the period from days to seconds
    total_mass = (mass1 + mass2)
    amplitude = ((((mass2**3) / (total_mass**2)) * (np.sin(inclination_angle)**3)
                  *((2 * pi * Gravitational_constant) / orbital_period))**(1/3))
    return (amplitude) #m/s

def calculate_radial_velocity(mass1, q_array, time_array, period, theta_array, inclination_angle_array):
    '''Here is a function to calculate the radial velocity of a star
    'q'is the mass2:mass1 ratio
    radial_velocity = Acos(wt+o),  where
    A=amplitude
    w=2*pi/period
    t=time
    o=theta
    period should be entered in days and the code will transform it to seconds,
    to get units of m/s'''
    time_array = time_array.reshape(len(time_array), 1, 1, 1)
    theta_array = theta_array.reshape(1, len(theta_array), 1, 1)
    q_array = q_array.reshape(1, 1, len(q_array),1)
    inclination_angle_array = inclination_angle_array.reshape(1, 1, 1, len(inclination_angle_array))
    amplitude = calculate_amplitude(mass1, q_array, inclination_angle_array, period) #m/s
    omega=2*pi/ period #period is in days
    #possible errors time arrange is to small, make it into bigger chuncks.
    radial_velocity = ((amplitude*np.cos((omega*time_array) + theta_array))) #units of m/s
    return radial_velocity #m/s

def plot_rv(mass1, q_array, time_array, period, theta_array, inclination_angle_array, number_of_observations):
    radial_velocity = calculate_radial_velocity(mass1, q_array, time_array, period, theta_array, inclination_angle_array)
    for x in range(0, len(theta_array)):
        plt.plot(time_array, radial_velocity[:, x], linestyle = ':')
        
    
    
    
    