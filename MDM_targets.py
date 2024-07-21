import astropy
import astropy.io.ascii as ascii
import astropy.io.fits
import numpy as np
from astropy.table import Table, join, Column
import pickle
from math import pi
import pandas as pd
import json
""" data_set1 has information about our target stars and we can use
    column 1 ( Kepler Input Catalog) to find our stars of interest.
    data_set2 is a second set of information about our targets and we will use 
    the Kepler Input Catalog to find our stars
"""
def readData():

    targets = Table.read('MDM_targets.csv', format = 'ascii.csv')
    targets.rename_column('ï»¿KIC', 'KIC')
    data_set1 = Table.read('ajab8a33t1_mrt.txt', format = 'ascii.mrt')
    data_set2 = Table.read('ajab8a33t2_mrt.txt', format = 'ascii.mrt')
    #This is a table that we will use for our research purposes
    McQuillan = Table.read('McQuillan.fit', format = 'fits')

    #In this line of code we are joining data set 1 data set 2 and our targets table
    table = join(data_set1, data_set2)
    target_table = join(targets, table)
    #This is the final table were we combine all of our data
    final_table = join(target_table, McQuillan, join_type = 'left', keys = 'KIC')
    #Now this are the tables of interes that we will use in our model
    return final_table

def rv_times():
    with open('full_rvs.pickle', 'rb') as f:
       docs, full_rvs = pickle.load(f) #pickle.load(f)
    KIC_times, KIC_rvs = full_rvs["KIC6425783"]
    return KIC_times, KIC_rvs

#Save full table to code
def writeData():
    full_table = readData()
    full_table.write('full_table.csv', format='ascii.csv')
  
#Read full table from the function above    
def readDataFromFile():
    full_table = Table.read('full_table.csv')
    return full_table

#Return target mass for an individual KIC number in SUN masses
def target_mass(KIC_number):
    full_table = readDataFromFile()
    target_mass = full_table['KIC', 'Mass']
    target_mass.add_index('KIC')
    mass = target_mass.loc[KIC_number]['Mass']
    return mass

#Return target period for an individual KIC number in days
def target_period(KIC_number):
    full_table = readDataFromFile()
    target_period = full_table['KIC', 'Prot']
    target_period.add_index('KIC')
    Prot = target_period.loc[KIC_number]['Prot']
    return Prot

'''
    This is a file with radial velocities, the are the actual velcites the number of 
    observation times depends on the size of the array
'''

#KIC numbers
def KIC_index():
    full_table = readDataFromFile()
    index =  full_table['KIC']
    return index
 
#This function returns a sin distribution of the inclination angles   
def sin_distn():
    inclination_angle_array = np.random.uniform(0, pi/2, 300)
    y = np.arccos((2 * inclination_angle_array)/pi)
    return y    

#A little bit of seth theory
def missingData():
    targets = readData().targets 
    table = readData().table
    A = set(targets['KIC'])
    B = set(table['KIC'])
    missing = A - (A & B)



