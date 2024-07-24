import astronomy_functions as af
import MDM_targets as md
import numpy as np
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

#Declare Variables and constants
SUN_MASS = const.M_sun # Value = 1.98840987e+30 kg

KIC_index = md.KIC_index()

def KIC_Model(KICnumber):
    #Declare paramters
    time_array, measured_rv = md.rv_times(KICnumber) 
    time_array = np.array(time_array) *u.day
    mass1 = md.target_mass(KICnumber) * SUN_MASS
    target_period = md.target_period(KICnumber)
    
    #Find minimum Chi-Square value and respective modeled parameters
    minimum = af.minimum(KICnumber)
    min_ = minimum
    b_f_theta = np.array([min_[0].value]) * u.radian
    b_f_q = np.array([min_[1]])
    b_f_i = np.array([min_[2]])
    min_chi_value = np.array([min_[3]])
    #print(b_f_theta,b_f_q,b_f_i)
    
    #calculate best fit curve radial velocity and plot it with measured rv
    b_f_rv = af.calculate_radial_velocity(mass1, b_f_q, time_array, target_period, b_f_theta, b_f_i)
    sns.set()
    ax = plt.axes()
    ax.set(xlabel = 'time (days)', ylabel = 'radial velocity (km/sec)',
           title = f"KIC_{KICnumber} observed & measured rvs")
    plt.plot(time_array, (b_f_rv[:,0,0,0].to(u.km/u.s)).value, label = 'best fit rv curve', linestyle = ':')
    plt.plot(time_array,measured_rv, label = 'measured rv', linestyle = '-')
    plt.legend()
    plt.show()
    
    #Hypothesis testing
    print('Experimental p value = ', min_chi_value)
    right_tail_value = stats.chi2.ppf(1-0.025, df = 2)#we use 2 degrees of freedom since we have 3 categories of modeled parameters, dof = n-1 = (theta, q, inc_angle) - 1
    print('Theoretical p value =', right_tail_value)
    if min_chi_value < right_tail_value:
        print('We have a binary system')
    elif min_chi_value > right_tail_value:
        print('We do not have a binary system')
    print("")

def Run_Model():
    for i in KIC_index:
        KIC_Model(i)