from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
from math import pi
import MDM_targets as md



#Declare parameter arrays
q_array = np.random.uniform(0,1, 20)
inclination_angle_array = md.sin_distn()
theta_array = np.random.uniform(0, 2 * pi, 10) * u.radian


Gravitational_constant = const.G # Value = 6.67 x 10^-11 m^3/(kg s^2)
def calculate_amplitude(mass1, q_array , inclination_angle_array, period): #m/s
    """ calculate the amplitude (semi-major axis) of the binary system. """
    mass2 = q_array * mass1 
    orbital_period = (period * u.day).to(u.second)
    total_mass = mass1 + mass2 # in kg
    amplitude = (
                    ( ((mass2**3) / (total_mass**2)) 
                    * (np.sin(inclination_angle_array)**3)
                  *((2 * pi * Gravitational_constant) / orbital_period))**(1/3)
                )
    return amplitude #units of m/s

def calculate_radial_velocity(mass1, q_array, time_array, period , theta_array, inclination_angle_array):#m/s
    """ Calculate the radial velocity of a star. """
    time_array = time_array.reshape(len(time_array), 1, 1, 1) 
    theta_array = theta_array.reshape(1, len(theta_array), 1, 1)
    q_array = q_array.reshape(1, 1, len(q_array),1)
    inclination_angle_array = inclination_angle_array.reshape(1, 1, 1, len(inclination_angle_array))
    
    amplitude = calculate_amplitude(mass1, q_array, inclination_angle_array, period) #m/s
    period = period * u.day
    omega=(2*pi) / (period.to(u.second)) #angular frequency in radians/s
    radial_velocity = ((amplitude*np.cos(((omega*time_array).to(u.radian, equivalencies = u.dimensionless_angles())) + theta_array))) #units of m/s
    return radial_velocity # in m/s


def calculate_velocity_with_noise(mass1, q_array, time_array, period, std, theta_array, inclination_angle_array):#km/s
    """ Calculate an array of radial velocities adding noise. """
    radial_velocity_array = calculate_radial_velocity(mass1, q_array, time_array, period,theta_array,inclination_angle_array) #m/s
    noise = np.random.normal(0, std, size = radial_velocity_array.shape)  * (u.km/u.second)
    velocity_with_noise = radial_velocity_array + noise
    return (velocity_with_noise.to(u.km/u.s)) # km/sec


SUN_MASS = const.M_sun # Value = 1.98840987e+30 kg
def generate_rv_array(KICnumber):
    """ Generate an array of radial velocity values for a given KIC target. """
    time_array, measured_rv = md.rv_times(KICnumber)
    time_array = np.array(time_array) *u.day
    measured_rv = np.array(measured_rv) * (u.km/u.s)
    std = 10 # standard deviation in km/s
    mass = md.target_mass(KICnumber) * SUN_MASS
    target_period = md.target_period(KICnumber) 
    rv_array = calculate_velocity_with_noise(mass, q_array, time_array, target_period, std, theta_array, inclination_angle_array) #m/sec
    return rv_array # in km/s


def Chi_square_distribution(KICnumber):
    """ Create a Chi-square value for each generated rv to compare it with the measured results. """
    radial_velocities = generate_rv_array(KICnumber)
    time_array,measured_rv = md.rv_times(KICnumber)
    measured_rv = measured_rv * (u.km/u.s) #measured rv
    measured_rv = measured_rv.reshape(len(measured_rv), 1, 1, 1)
    chi_square = ((radial_velocities - measured_rv)**2) / ((10* u.km/u.s)**2)
    chi_values_array = np.sum(chi_square, axis = 0)
    return chi_values_array

              
def minimum(KICnumber):
    """ Find the minimum Chi-square vale to find the parameters for the best-fit radial velocity curve. """
    chi_values = Chi_square_distribution(KICnumber)
    min_value = chi_values.min() #returns the minimum value in the multidimensional array s
    max_value = chi_values.max()
    min_index = np.argwhere(chi_values == min_value)
    
    theta_value = theta_array[min_index[0][0]]
    q_value = q_array[min_index[0][1]]
    i_value = inclination_angle_array[min_index[0][2]]
    
    print(f"KIC target: KIC_{KICnumber}")
    print("max value is ", max_value)
    print("min value is ", min_value)
    print("theta value is ", theta_value)
    print("q value is ", q_value)
    print("inclination angle value is ", i_value)
    
    return [theta_value,q_value,i_value,min_value]
      

# fig.savefig('my_figure.png') saves figures might use it to save plots to specific locations
# from IPython.display import Image
# Image('my_figure.png') displays a saved image or plot as defined above

    