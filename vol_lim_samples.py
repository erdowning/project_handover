import astropy
from astropy.table import Table
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as U

import numpy as np
import matplotlib.pyplot as plt
import os 

from smith_kcorr import GAMA_KCorrection


cosmo = FlatLambdaCDM(H0=100, Om0=0.3)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def absmag(appmag, z, k, e):
    return appmag - 5*np.log10(cosmo.luminosity_distance(np.array(z)).to_value()) - 25 - k - e

def appmag(absmag, z, k, e):
    return absmag + 5*np.log10(cosmo.luminosity_distance(np.array(z)).to_value()) + 25 + k + e

def vol_lim_sample_2(lim0, lim1, cat, apply_absmask=True): #lim0=faintlim, lim1=brightlim
    kcorr_r_N = GAMA_KCorrection(band='R', file='jmext', photsys='N')
    kcorr_r_S = GAMA_KCorrection(band='R', file='jmext', photsys='S')
    
    N_mask = (cat['PHOTSYS'] == 'N')
    S_mask = (cat['PHOTSYS'] == 'S')
    cat['KMIN'] = np.zeros(len(cat))
    cat['KMAX'] = np.zeros(len(cat))
    cat['KMIN'][N_mask] = kcorr_r_N.k(cat['Z'][N_mask], np.ones(len(cat['Z'][N_mask]))*0.39861392974853516) # k-correction for galaxy's redshift but bluest colour
    cat['KMAX'][N_mask] = kcorr_r_N.k(cat['Z'][N_mask], np.ones(len(cat['Z'][N_mask]))*0.9832057952880859) # k-correction for galaxy's redshift but reddest colour (for conservative faint  limit)
    cat['KMIN'][S_mask] = kcorr_r_S.k(cat['Z'][S_mask], np.ones(len(cat['Z'][S_mask]))*0.4014892578125) #NEEDS TO BER CHANGED WITH K CORRECTIONS
    cat['KMAX'][S_mask] = kcorr_r_S.k(cat['Z'][S_mask], np.ones(len(cat['Z'][S_mask]))*0.9937572479248047)                                                                                                          
    
    cat['ABS_BRIGHTLIM'] =  absmag(15, np.array(cat['Z']), np.array(cat['KMIN']), np.array(cat['EQ_ALL_0P1']))
    cat['ABS_FAINTLIM'] = np.zeros(len(cat))
    cat['ABS_FAINTLIM'][N_mask] =  absmag(19.539993, np.array(cat['Z'][N_mask]), np.array(cat['KMAX'][N_mask]), np.array(cat['EQ_ALL_0P1'][N_mask])) #cat['FAINT_RLIM']
    cat['ABS_FAINTLIM'][S_mask] =  absmag(19.5, np.array(cat['Z'][S_mask]), np.array(cat['KMAX'][S_mask]), np.array(cat['EQ_ALL_0P1'][S_mask]))
    
    if lim1 != None and lim1 != 'None':
        volmask = ((lim0 <= cat['ABS_FAINTLIM']) & (lim1 >= cat['ABS_BRIGHTLIM']))

        if apply_absmask:
            try:
                absmask = ((cat['ABSMAG_R'] >= lim1) & (cat['ABSMAG_R'] < lim0)) #sv3

            except KeyError:
                absmask = ((cat['ABSMAG_RP1'] >= lim1) & (cat['ABSMAG_RP1'] < lim0)) #da0.2
            mask = (absmask & volmask)
        else:
            mask = volmask
            
    else:
        print(lim1)
        volmask = (lim0 <= cat['ABS_FAINTLIM']) & (cat['Z'] > 0.05) # magnitude threshold samples
        if apply_absmask:
            try:
                absmask = (cat['ABSMAG_R'] < lim0) #sv3

            except KeyError:
                absmask = (cat['ABSMAG_RP1'] < lim0) #da0.2
            mask = (absmask & volmask)
        else:
            mask = volmask
        
    return mask