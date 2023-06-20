#This code is used to determine the influence of
#binary stars on the CMD
#The idea is to take an isochrone which fits the
#observational data well and use it to sample 
#binary stars with different qs. The result is 
#our understanding of the morphology of binary
#stars on the CMD which can be used to sample
#binary stars with the fiducial isochrones

import numpy as np
import pandas as pd
import os
from vorbin_chi2.utiles import check_file
from scipy.interpolate import interp1d

def read_iso(iso_path,age):
    check_file(iso_path)
    iso = open(iso_path, 'r')
    len_file = len(iso.readlines())
    iso.seek(0)
    AGE_wanted = np.float(age)
    retval = []
    if len_file < 10:
        raise Exception("Empty file")
    else:   
        #skip header
        iso.readline()
        NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
        while int(AGE[:-1]) != AGE_wanted:
            #skip header
            iso.readline()
            for i in range(int(NPTS[1:])):
                #skiplines
                iso.readline()
            #skiplines
            iso.readline()
            iso.readline()
            iso.readline()
            NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
        #skip header
        iso.readline()
        for i in range(int(NPTS[1:])):
            EEP,MASS,LogG,LogTeff,LogL,_,_,_,V,VI,F606W,F606WF814W = iso.readline().split()
            retval.append([EEP,MASS,LogG,LogTeff,LogL,V,VI,F606W,F606WF814W])
    dp = pd.DataFrame(data=retval,columns=['EEP','MASS','LogG','LogTeff','LogL','V','VI','F606W','F606WF814W'])
    return dp

def find_MSTO(iso):
    vi = iso['F606WF814W'].values
    v = iso['F606W'].values
    #define the bluest point to be MSTO
    vi_min = np.inf
    for i in range(len(iso)):
        if vi[i] < vi_min:
            vi_min = vi[i]
            MSTO_idx = i
    return MSTO_idx, v[MSTO_idx], vi[MSTO_idx]

def find_SGB(iso, MSTO_idx, MSTO_V, MSTO_VI):
    vi = iso['F606WF814W'].values
    v = iso['F606W'].values
    #define the Sub Giant Branch to be 0.05 mag redder from MSTO
    SGB_vi = MSTO_VI + 0.05
    SGB_v = np.interp(SGB_vi, vi[MSTO_idx:], v[MSTO_idx:])
    return SGB_v, SGB_vi

def define_bound(iso, mag_cut=2):
    MSTO_idx, MSTO_V, MSTO_VI = find_MSTO(iso)
    SGB_V, SGB_VI = find_SGB(iso, MSTO_idx, MSTO_V, MSTO_VI)
    return SGB_V-mag_cut, SGB_V+mag_cut

def define_func(iso):
    vi = iso['F606WF814W'].values
    v = iso['F606W'].values
    m = iso['MASS'].values
    m_func = interp1d(v,m,kind='‘cubic’')
    v_func = interp1d(m,v,kind='‘cubic’')
    vi_func = interp1d(m,vi,kind='‘cubic’')
    return m_func, v_func, vi_func

def __init__(iso_path, age, mag_cut=2):
    iso = read_iso(iso_path, age)
    v_min, v_max = define_bound(iso, mag_cut=mag_cut)
    m_func, v_func, vi_func = define_func(iso)