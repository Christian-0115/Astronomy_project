import astropy
import astropy.io.ascii as ascii
import astropy.io.fits
import numpy as np
from astropy.table import Table, join, Column
import pickle
from math import pi

def readData():
    """ Read and merge the necessary data sets into a final table. """
    targets = Table.read('MDM_targets.csv', format = 'ascii.csv')
    targets.rename_column('ï»¿KIC', 'KIC')
    data_set1 = Table.read('ajab8a33t1_mrt.txt', format = 'ascii.mrt')
    data_set2 = Table.read('ajab8a33t2_mrt.txt', format = 'ascii.mrt')
    McQuillan = Table.read('McQuillan.fit', format = 'fits')

    # Join the data sets 
    table = join(data_set1, data_set2)
    target_table = join(targets, table)
    final_table = join(target_table, McQuillan, join_type = 'left', keys = 'KIC')
    
    return final_table

def rv_times(KIC_number):
    """ Get radial velocity values and observation times for a given KIC number. """
    with open('full_rvs.pickle', 'rb') as f:
       docs, full_rvs = pickle.load(f) 
    KIC_times, KIC_rvs = full_rvs[f"KIC{KIC_number}"]
    
    return KIC_times, KIC_rvs

def writeData():
    """ Save the final table into a CSV file"""
    full_table = readData()
    full_table.write('full_table.csv', format='ascii.csv') 

def readDataFromFile():
    """ Read CSV file full table. """
    return Table.read('full_table.csv')
   
def target_mass(KIC_number):
    """ Return target mass for an individual KIC number in SUN masses. """
    full_table = readDataFromFile()
    target_mass = full_table['KIC', 'Mass']
    target_mass.add_index('KIC')
    mass = target_mass.loc[KIC_number]['Mass']
    return mass

def target_period(KIC_number):
    """ Return target period for an individual KIC number in days. """
    full_table = readDataFromFile()
    target_period = full_table['KIC', 'Prot']
    target_period.add_index('KIC')
    Prot = target_period.loc[KIC_number]['Prot']
    
    return Prot

#KIC numbers
def KIC_index():
    """ Return the KIC index of all our targets from full table. """
    full_table = readDataFromFile()
    index =  full_table['KIC']
    return index
   
def sin_distn():
    """ Generate a sine distribution of the inclination angles. """
    inclination_angle_array = np.random.uniform(0, pi/2, 30)
    y = np.arccos((2 * inclination_angle_array)/pi)
    return y    

def missingData():
    """ Identify missing data in the target table. """
    targets = readData().targets 
    table = readData().table
    A = set(targets['KIC'])
    B = set(table['KIC'])
    missing = A - (A & B)
    return missing



