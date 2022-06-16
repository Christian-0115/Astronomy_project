# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import matplotlib.pyplot as plt
import numpy as np
from math import pi
import seaborn as sns

EARTH_MASS = 6 * (10**24)
SUN_MASS = 2 * (10**30)

def calculate_amplitude(mass1, q , inclination_angle, period):
    """ Here is another function to calculate the amplitude of a binary system
    in a binary system.
    The inputs consist of:
        mass of star1 (bigger object)
        'q' which is the mass2:mass1 ratio 
        angle of inclination
        and period of orbit in days, and the code will change it to seconds
        the amplitude output has units of km/sec
    """
    mass2 = q * mass1
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
    return radial_velocity


def calculate_velocity_with_noise(mass1, q_array, time_array, period, std, theta_array,
                                  inclination_angle_array, number_of_observations):
    #our code is correct up to this point it'll be useful to check it from here on
    """ Here is a function to calculate an array of velocities with noise
        period should be entered in days to get units of m/s
        standard deviation needs to have units of m/sec
        'q'is the mass2:mass1 ratio
    """
    radial_velocity_array = calculate_radial_velocity(mass1, q_array, time_array, period,
                                                      theta_array,inclination_angle_array) #m/s
    noise = np.random.normal(0, std, number_of_observations) #km/sec
    noise_array = noise.reshape(len(noise), 1, 1, 1)
    velocity_with_noise = radial_velocity_array + noise_array
    return velocity_with_noise #m/sec

def plot_velocity_vs_time(mass1, q, period, time_array, theta_array, inclination_angle, std, 
                          number_of_observations):
    """ Here we are creating a function to graph the radial velocity
        of a binary system, when we are given:
            'q'is the mass2:mass1 ratio
            set the time array as np.linspace (0, 2*period, #of data points)
            set the theta array as np.random.uniform(0, 2*pi, # of angles desired)
            mass1, mass2, period, std, theta, and the inclination angle.
        The function returns velocity vs time over 2 orbital periods.
    """
    sns.set()
    radial_velocity_with_noise = calculate_velocity_with_noise(mass1, q, time_array, 
                                                          period, std, theta_array, inclination_angle, 
                                                          number_of_observations) #m/sec
    ax = plt.axes()
    ax.set(xlabel = 'time (days)', ylabel = 'radial velocity (km/sec)',
           title = 'Motion of star 1 around star 2')
    for x in range(0, len(theta_array)):
        plt.plot(time_array, (radial_velocity_with_noise[x, :] / 1000), linestyle = ':') #km/s
        plt.errorbar(time_array, (radial_velocity_with_noise[x, :] / 1000), yerr = std, fmt = '.',
                     ecolor = 'blue');

def Chi_square_distribution (mass1, q, period, time_array, theta_array, inclination_angle,
                             std, number_of_observations):
    """ In this function we will be running a Chi Square statistical test to see if
        if any of our systems fall within the range of a binary star distribution
    """
    
    
    
    return

           
def plot_velocity_vs_time_with_histogram(mass1, q, time_array, period, theta_array, 
                                         inclination_angle, std, number_of_observations):
    ''' period has to be in seconds
        units of radial velocity are in km/sec
        'q'is the mass2:mass1 ratio
    ''' 
    radial_velocity_with_noise = calculate_velocity_with_noise(mass1, q, time_array, 
                                                          period, std, theta_array, inclination_angle, 
                                                          number_of_observations)    
    fig = plt.figure()
    ax1 = fig.add_axes()
    ax2 = fig.add_axes()
    ax1 = plt.plot(time_array, radial_velocity_with_noise, linestyle = ':', color = 'black')
    ax2 = plt.hist(radial_velocity_with_noise)
    return ax1, ax2

def test_function():
    print("I want to test this!")
    print("Here is another test")

def another_test():
    print("I am doing another test.")

 



    
    