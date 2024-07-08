import astropy
import astropy.io.ascii as ascii
import astropy.io.fits
import numpy as np
from astropy.table import Table, join, Column
import pickle
from math import pi
import pandas as pd

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

#table info ^^^
''' table info ^^
 name    dtype      unit                      description                      class     n_bad 
-------- ------- ------------ ---------------------------------------------- ------------ ------
     KIC   int32                             Kepler Input Catalog identifier       Column      0
    gmag float64          mag                           g band magnitude (1)       Column      0
  e_gmag float64          mag                            Uncertainty in gmag       Column      0
   Ksmag float64          mag                          Ks band magnitude (2)       Column      0
 e_Ksmag float64          mag                           Uncertainty in Ksmag       Column      0
     Par float64          mas                              Gaia DR2 Parallax       Column      0
   e_Par float64          mas                             Uncertainty in Par       Column      0
  [Fe/H] float64     dex(---)                  Spectroscopic metallicity (3) MaskedColumn 116841
e_[Fe/H] float64     dex(---)                          Uncertainty in [Fe/H] MaskedColumn 116841
    RUWE float64                         Re-normalized unit-weight error (4)       Column      0
   Ncomp   int32                       Number, Gaia DR2 companions within 4"       Column      0
  KsCorr   str14              Potential corrections compared to 2MASS Ks (5) MaskedColumn 152906
   State    str5                                      Evolutionary State (6) MaskedColumn 171509
   Mstar float64         Msun    Isochrone derived stellar mass, solar units       Column      0
 E_Mstar float64         Msun                           Upper error on Mstar       Column      0
 e_Mstar float64         Msun                           Lower error on Mstar       Column      0
    Teff float64            K        Isochrone derived effective temperature       Column      0
  E_Teff float64            K                            Upper error on Teff       Column      0
  e_Teff float64            K                            Lower error on Teff       Column      0
    logg float64 dex(cm / s2)         Isochrone derived surface gravity, log       Column      0
  E_logg float64 dex(cm / s2)                            Upper error on logg       Column      0
  e_logg float64 dex(cm / s2)                            Lower error on logg       Column      0
     FeH float64     dex(---)          Isochrone derived surface metallicity       Column      0
   E_FeH float64     dex(---)                             Upper error on FeH       Column      0
   e_FeH float64     dex(---)                             Lower error on FeH       Column      0
   Rstar float64         Rsun  Isochrone derived stellar radius, solar units       Column      0
 E_Rstar float64         Rsun                           Upper error on Rstar       Column      0
 e_Rstar float64         Rsun                           Lower error on Rstar       Column      0
     rho float64     dex(---)    Isochrone derived density, solar units, log       Column      0
   E_rho float64     dex(---)                             Upper error on rho       Column      0
   e_rho float64     dex(---)                             Lower error on rho       Column      0
   Lstar float64    dex(Lsun) Isochrone derived luminosity, solar units, log       Column      0
 E_Lstar float64    dex(Lsun)                           Upper error on Lstar       Column      0
 e_Lstar float64    dex(Lsun)                           Lower error on Lstar       Column      0
     Age float64          Gyr                     Isochrone derived age, Gyr       Column      0
   f_Age    str1                                           [* ] Age flag (1) MaskedColumn 159720
   E_Age float64          Gyr                             Upper error on Age       Column      0
   e_Age float64          Gyr                             Lower error on Age       Column      0
    Dist float64           pc                 Isochrone derived distance, pc       Column      0
  E_Dist float64           pc                            Upper error on Dist       Column      0
  e_Dist float64           pc                            Lower error on Dist       Column      0
   Avmag float64          mag            Isochrone derived V-band extinction       Column      0
     GOF float64                         combined likelihood goodness-of-fit       Column      0
    TAMS float64          Gyr              Terminal age of the main sequence       Column      0
'''

'''
    This file contains a tuple consisting of these instructions as the first \nelement,
    and a dictionary mapping the KIC name of the targets to a 2-tuple\nconaining the HJD 
    times of the observation and the measured RVs.
'''
def rv_times():
    with open('full_rvs.pickle', 'rb') as f:
        full_rvs = pickle.load(f)
    return full_rvs

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
def observation_times(KIC_number):
    observation_times = rv_times()
    return observation_times

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



