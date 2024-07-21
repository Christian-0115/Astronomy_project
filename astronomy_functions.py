from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from math import pi
import seaborn as sns
import scipy.stats as stats
import MDM_targets as md
import pandas as pd

#Declare constants 
EARTH_MASS = const.M_earth # Value = 5.97216787e+24 kg
SUN_MASS = const.M_sun # Value = 1.98840987e+30 kg


#Declare parameter arrays
q_array = np.random.uniform(0,1, 200)
inclination_angle_array = md.sin_distn()
theta_array = np.random.uniform(0, 2 * pi, 100) * u.radian

time_array, measured_rv = md.rv_times()
time_array = np.array(time_array) *u.day
measured_rv = np.array(measured_rv) * (u.km/u.s)

#This is a function to calculate the semi amplitude of RV - note need to researhch what this function actually is
def calculate_amplitude(mass1, q_array , inclination_angle_array, period): #m/s
    """ 
    Here is a function to calculate an array of amplitudes for a binary system, 
    the size of the array depends on the size of the input arrays
    The inputs consist of:
        mass of star1 (bigger object) in kg.
        'q' which is the mass2:mass1 ratio 
        angle of inclination
        and period of orbit in days (the code will change it to seconds)
        the amplitude output has units of m/sec
    """
    mass1 = (mass1)# mass1 is in kg
    mass2 = q_array * (mass1) # mass 2 is in kg
    Gravitational_constant = const.G # Value = 6.67 x 10^-11 m^3/(kg s^2)
    orbital_period = period * u.day
    orbital_period = orbital_period.to(u.second)#change the period from days to seconds
    total_mass = (mass1 + mass2) # in kg
    amplitude = (
                    ( ((mass2**3) / (total_mass**2)) 
                    * (np.sin(inclination_angle_array)**3)
                  *((2 * pi * Gravitational_constant) / orbital_period))**(1/3)
                )
    return amplitude #units of m/s

'''Here is a function to calculate the radial velocity of a star
    q_array is the mass2:mass1 ratio
    radial_velocity = Acos(wt+o),  where
    A=amplitude in m/s
    w=2*pi/period in days
    t=time_array
    o=theta_array
    period should be entered in days and the code will transform it to seconds,
    to get units of m/s
'''
def calculate_radial_velocity(mass1, q_array, time_array, period , theta_array, inclination_angle_array):#m/s
    #Specify inputs as quantity objects
    time_array = time_array.reshape(len(time_array), 1, 1, 1) # time array 
    theta_array = theta_array.reshape(1, len(theta_array), 1, 1) #Phase angle array
    q_array = q_array.reshape(1, 1, len(q_array),1)
    inclination_angle_array = inclination_angle_array.reshape(1, 1, 1, len(inclination_angle_array))
    amplitude = calculate_amplitude(mass1, q_array, inclination_angle_array, period) #m/s
    period = period * u.day
    omega=(2*pi) / (period.to(u.second)) #period is in seconds
    radial_velocity = ((amplitude*np.cos(((omega*time_array).to(u.radian, equivalencies = u.dimensionless_angles())) + theta_array))) #units of m/s
    return radial_velocity # in m/s

""" Here is a function to calculate an array of velocities with noise
    period should be entered in days to get units of m/s
    standard deviation needs to have units of m/sec
    'q'is the mass2:mass1 ratio
"""
def calculate_velocity_with_noise(mass1, q_array, time_array, period, std, theta_array,
                                  inclination_angle_array):#km/s
    radial_velocity_array = calculate_radial_velocity(mass1, q_array, time_array, period,theta_array,inclination_angle_array) #m/s
    noise = np.random.normal(0, std, size = radial_velocity_array.shape)  * (u.km/u.second)
    velocity_with_noise = radial_velocity_array + noise
    return (velocity_with_noise.to(u.km/u.s)) #km/sec


'''
    This function generates an array of radial velocity values using known quantities about the stars
    and some of the modeled quantities. This array is pretty big and we will use it in a chi-square
    test to find the best fit radial velocity curve for a given KIC target

9653110
9655045
9710336
9964938
10153521
11819949
12736892

'''
def generate_rv_array(KICnumber):
    std = 10 #our standard deviation will be 10 km/s
#    time_array = KIC_6425783_time_of_observation*u.day # this gives us days of observation over a certain number of periods
    mass = md.target_mass(KICnumber) * SUN_MASS
    target_period = md.target_period(KICnumber) #* u.day
    inclination_angle_array = md.sin_distn()
    rv_array = calculate_velocity_with_noise(mass, q_array, time_array, target_period, std, theta_array, inclination_angle_array) #m/sec
    return rv_array # output in km/s


""" In this function we will be running a Chi Square statistical test to see if
    if any of our target stars fall within the range of a binary star distribution.
    In other words, this function will generate a chi-square value for each radial velocity curve
    in the array, then we will use the function after to find the minimum Chi-squre value
    which corresponds to the best-fit curve
    

9653110
9655045
9710336
9964938
10153521
11819949
12736892

""" #need to work out the units
def Chi_square_distribution(KICnumber):
    radial_velocities = generate_rv_array(KICnumber)
    time_array,measured_rv = md.rv_times()
    measured_rv = measured_rv * (u.km/u.s) #measured rv
    measured_rv = measured_rv.reshape(len(measured_rv), 1, 1, 1)
    chi_square = ((radial_velocities - measured_rv)**2) / ((10* u.km/u.s)**2)
    chi_values_array = np.sum(chi_square, axis = 0)
    return chi_values_array

      
'''
    This function is to find the minimum Chi-square value which corresponds tp the 
    best-fit radial velocity curve for a given KIC target with the given parameters
'''          
def minimum(KICnumber):
    s = Chi_square_distribution(KICnumber)
    min_value = s.min() #returns the minimum value in the multidimensional array s
    max_value = s.max()
    print("max value is ", max_value)
    print("min value is ", min_value)
    min_index = np.argwhere(s == min_value)
    # print("The index of the minimum value is ", min_index)
    theta_value = theta_array[min_index[0][0]]#0,0 specifies the value of theta and all rows
    q_value = q_array[min_index[0][1]]#0,1 specifies the value of q and all rows
    i_value = inclination_angle_array[min_index[0][2]]#0,1 specifies the value of q and all rows
    print("theta value is ", theta_value)
    print("q value is ", q_value)
    print("inclination angle value is ", i_value)
    best_fit_parameters = [theta_value,q_value,i_value,min_value]
    return best_fit_parameters#This returns an array for best-fit parameters of our model
      
def hist():
    df = theta_array.value
    plt.hist(df, bins = 500, color = 'slategrey', normed = 1)
    plt.title("Phase Angle Probability Distribution")
    

    

    
#np.argwhere(s) returns the shape of the matrix
# fig.savefig('my_figure.png') saves figures might use it to save plots to specific locations
# from IPython.display import Image
# Image('my_figure.png') displays a saved image or plot as defined above

# Have a graph that will output a single line for any single combination of inputs minus time time will always be constant
# if you have an array you have to go to speed.to()value and it will take the array out of the function
    