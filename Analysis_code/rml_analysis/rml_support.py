#This file includes all the support functions for RML analysis

import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import os
from typing import Union
import numpy.typing as npt
FARRAY_1D = npt.NDArray[np.float64]
def check_file(path):
    if os.path.exists(path) == False:
        raise Exception("Cannot find inputfile at {}".format(path))

def read_candidates(star_path):
    dp = pd.read_csv(star_path)
    N_stars = len(dp)
    Names = dp['Names'].values
    Mass_star = dp['M/Ms'].values
    Mass_star_p_err = dp['+dM/Ms'].values
    Mass_star_m_err = dp['-dM/Ms'].values
    Lumi_star = dp['L/Ls'].values
    Lumi_star_p_err = dp['+dL/Ls'].values
    Lumi_star_m_err = dp['-dL/Ls'].values
    Rad_star = dp['R/Rs'].values
    Rad_star_p_err = dp['+dR/Rs'].values
    Rad_star_m_err = dp['-dR/Rs'].values
    return Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err

def calculate_chi2( df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, Consider_Lumi=True):
    N_stars = len(Mass_star)
    #define the matrix value
    age_num = len(df_iso.index.get_level_values(0).unique())
    chi2_data = np.zeros((age_num, N_stars))
    if Consider_Lumi == True:
        for i in range(N_stars):
            Mass = df_iso.loc[(slice(None),'Mass'),:].copy().values
            Lumi = df_iso.loc[(slice(None),'Lumi'),:].copy().values
            Rad = df_iso.loc[(slice(None),'Rad'),:].copy().values
            #calculate difference
            Mass -= Mass_star[i]
            Lumi -= Lumi_star[i]
            Rad -= Rad_star[i]
            #find whether the difference is greater than 0 or not
            Mass_TF = Mass > 0
            Lumi_TF = Lumi > 0
            Rad_TF = Rad > 0
            #for each star, evaluate whether the difference is positive or not. Then divide by the corresponding uncertainty
            Mass_chi2 = np.divide(Mass,np.where(Mass_TF,Mass_star_p_err[i],Mass_star_m_err[i]))
            Lumi_chi2 = np.divide(Lumi,np.where(Lumi_TF,Lumi_star_p_err[i],Lumi_star_m_err[i]))
            Rad_chi2 = np.divide(Rad,np.where(Rad_TF,Rad_star_p_err[i],Rad_star_m_err[i]))
            #square the result
            Mass_chi2 = np.square(Mass_chi2)
            Lumi_chi2 = np.square(Lumi_chi2)
            Rad_chi2 = np.square(Rad_chi2)
            chi2 = Mass_chi2 + Lumi_chi2 + Rad_chi2
            indi_fit = np.min(chi2,axis=1)
            chi2_data[:,i] = indi_fit
    else:
        for i in range(N_stars):
            Mass = df_iso.loc[(slice(None),'Mass'),:].copy().values
            Rad = df_iso.loc[(slice(None),'Rad'),:].copy().values
            #calculate difference
            Mass -= Mass_star[i]
            Rad -= Rad_star[i]
            #find whether the difference is greater than 0 or not
            Mass_TF = Mass > 0
            Rad_TF = Rad > 0
            #for each star, evaluate whether the difference is positive or not. Then divide by the corresponding uncertainty
            Mass_chi2 = np.divide(Mass,np.where(Mass_TF,Mass_star_p_err[i],Mass_star_m_err[i]))
            Rad_chi2 = np.divide(Rad,np.where(Rad_TF,Rad_star_p_err[i],Rad_star_m_err[i]))
            #square the result
            Mass_chi2 = np.square(Mass_chi2)
            Rad_chi2 = np.square(Rad_chi2)
            chi2 = Mass_chi2 + Rad_chi2
            indi_fit = np.min(chi2,axis=1)
            chi2_data[:,i] = indi_fit
    return chi2_data

def writeout(dp,wrt_path):
    dp.to_csv(wrt_path,index=False,mode='a',header=not os.path.exists(wrt_path))

def read_iso(iso_path, max_eep=1500 ,age_list: Union[bool, FARRAY_1D] = np.linspace(8000,16000,81).astype(int),param_list=np.array(['Mass', 'Lumi', 'Rad']), iso_D = False):
    check_file(iso_path)
    mc_num = int(iso_path[-5:])
    iso = open(iso_path, 'r')
    len_file = len(iso.readlines())
    iso.seek(0)
    len_idx = 1
    if len_file < 10:
        iso.close()
        raise Exception("Empty file")
    else:
        if iso_D == False:
            if age_list == 'Unknown':
                age_list = []
                while len_idx <= len_file:
                    #look for how many ages it contains
                    iso.readline()
                    RAW_NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
                    npts = int(RAW_NPTS[1:])
                    age_list.append(int(AGE[:-1]))
                    iso.readline()
                    len_idx += 3
                    for i in range(npts):
                        iso.readline()
                    len_idx += npts
                    iso.readline()
                    iso.readline()
                    len_idx += 2
                #return to default
                iso.seek(0)
                len_idx = 1
                age_list = np.array(age_list)
            #generate Mutiindex dataframe
            num_age = len(age_list)
            num_param = len(param_list)
            arrays = [np.repeat(age_list,num_param),np.repeat(param_list.reshape(1,num_param),num_age, axis=0).reshape(num_age*num_param)]
            data = np.full((num_age*num_param, max_eep),np.inf)
            data_idx = 0
            while len_idx <= len_file:
                #skip header
                iso.readline()
                RAW_NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
                npts = int(RAW_NPTS[1:])
                # Mass = np.zeros(npts)
                # Lumi = np.zeros(npts)
                # Rad = np.zeros(npts)
                #skip header
                iso.readline()
                len_idx += 3
                for i in range(npts):
                    EEP,MMs,LogG,LogTeff,LogLLs,LogRRs = iso.readline().split()
                    # Mass[i] = float(MMs)
                    # Lumi[i] = float(LogLLs)
                    # Rad[i] = float(LogRRs)
                    data[data_idx,i] = float(MMs)
                    data[data_idx+1,i] = 10**float(LogLLs)
                    data[data_idx+2,i] = 10**float(LogRRs)
                data_idx += 3
                len_idx += npts
                iso.readline()
                iso.readline()
                len_idx += 2
            iso.close()
        else:
            iso.close()
            #generate Mutiindex dataframe
            num_age = len(age_list)
            num_param = len(param_list)
            arrays = [np.repeat(age_list,num_param),np.repeat(param_list.reshape(1,num_param),num_age, axis=0).reshape(num_age*num_param)]
            data = np.full((num_age*num_param, max_eep),np.inf)
            #read in isochrone file
            df = pd.read_csv(iso_path)
            df['Age'] = (df['Age']/10**6).astype(int)
            for i, age in enumerate(age_list):
                data[3*i,df[df['Age'] == age]['EEP']] = df[df['Age'] == age]['MMs']
                data[3*i+1,df[df['Age'] == age]['EEP']] = 10**(df[df['Age'] == age]['LogLLs'])
                data[3*i+2,df[df['Age'] == age]['EEP']] = 10**(df[df['Age'] == age]['LogRRs'])
        df = pd.DataFrame(data,index=arrays)
    return mc_num, df, age_list

def resample_stars( df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err):
    #given the information about an isochrone, and the observed binaries, we redraw stars from the theoretical isochrone with know uncertainty.
    num_stars = len(Mass_star)
    iso_Mass = df_iso.loc['Mass'].values
    TF = iso_Mass < 1e10
    #Remove inf
    iso_Mass = iso_Mass[TF]
    iso_Lumi = df_iso.loc['Lumi'].values[TF]
    iso_Rad = df_iso.loc['Rad'].values[TF]
    NPT = len(iso_Mass)
    pick_idx = np.random.choice(range(int(NPT/4),NPT - int(NPT/4)), size=num_stars, replace=False)
    Mass_alter = iso_Mass[pick_idx]
    Lumi_alter = iso_Lumi[pick_idx]
    Rad_alter = iso_Rad[pick_idx]
    sample_error= pd.DataFrame(columns = ['+dM/Ms','-dM/Ms','+dL/Ls','-dL/Ls','+dR/Rs','-dR/Rs'], data = np.abs(np.random.normal(scale=np.array([Mass_star_p_err, Mass_star_m_err, Lumi_star_p_err, Lumi_star_m_err, Rad_star_p_err, Rad_star_m_err]).T)))
    sample_error[['-dM/Ms','-dL/Ls','-dR/Rs']] = -sample_error[['-dM/Ms','-dL/Ls','-dR/Rs']]
    #return mass lumi and rad for resampled stars
    return Mass_alter + np.array([np.random.choice(sample_error[['+dM/Ms','-dM/Ms']].values[i]) for i in range(num_stars)]), Lumi_alter + np.array([np.random.choice(sample_error[['+dL/Ls','-dL/Ls']].values[i]) for i in range(num_stars)]), Rad_alter + np.array([np.random.choice(sample_error[['+dR/Rs','-dR/Rs']].values[i]) for i in range(num_stars)])

def resample_stars(df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err,EEPs, resample_num,Uniform_Resample=True):
    #given the information about an isochrone, and the observed binaries, we redraw stars from the theoretical isochrone with know uncertainty.
    #If Uniform_Resample == True, we draw stars uniformly from the center 50% of the isochrones
    #Otherwise, we choose resample stars from +/- 3 eeps from the given eep point
    num_stars = len(Mass_star)
    iso_Mass = df_iso.loc['Mass'].values
    iso_Lumi = df_iso.loc['Lumi'].values
    iso_Rad = df_iso.loc['Rad'].values
    if Uniform_Resample == True:
        iso_TF = iso_Mass < np.inf
        iso_Mass = iso_Mass[iso_TF]
        iso_Lumi = iso_Lumi[iso_TF]
        iso_Rad = iso_Rad[iso_TF]
        pick_idx = np.random.choice(range(int(len(iso_Mass)/4),len(iso_Mass) - int(len(iso_Mass)/4)), size=(num_stars, resample_num), replace=True)
    else:
        pick_idx = np.array([np.random.choice(range(EEPs[i-3],EEPs[i+3]), size=(num_stars, resample_num), replace=True) for i in range(num_stars)])
    Mass_alter = iso_Mass[pick_idx]
    Lumi_alter = iso_Lumi[pick_idx]
    Rad_alter = iso_Rad[pick_idx]
    M_p_re = np.abs(np.random.normal(scale=np.repeat(Mass_star_p_err.reshape(num_stars,1),resample_num,axis=1)))
    M_m_re = - np.abs(np.random.normal(scale=np.repeat(Mass_star_m_err.reshape(num_stars,1),resample_num,axis=1)))
    M_TF = np.random.choice(a=[1, 0], size=(num_stars, resample_num))
    M_err = np.multiply(M_p_re,M_TF) + np.multiply(M_m_re, (M_TF - 1)*-1)
    R_p_re = np.abs(np.random.normal(scale=np.repeat(Rad_star_p_err.reshape(num_stars,1),resample_num,axis=1)))
    R_m_re = - np.abs(np.random.normal(scale=np.repeat(Rad_star_m_err.reshape(num_stars,1),resample_num,axis=1)))
    R_TF = np.random.choice(a=[1, 0], size=(num_stars, resample_num))
    R_err = np.multiply(R_p_re,R_TF) + np.multiply(R_m_re, (R_TF - 1)*-1)
    L_p_re = np.abs(np.random.normal(scale=np.repeat(Lumi_star_p_err.reshape(num_stars,1),resample_num,axis=1)))
    L_m_re = - np.abs(np.random.normal(scale=np.repeat(Lumi_star_m_err.reshape(num_stars,1),resample_num,axis=1)))
    L_TF = np.random.choice(a=[1, 0], size=(num_stars, resample_num))
    L_err = np.multiply(L_p_re,L_TF) + np.multiply(L_m_re, (L_TF - 1)*-1)
    #return mass lumi and rad for resampled stars
    return (Mass_alter + M_err).T, (Lumi_alter + L_err).T, (Rad_alter + R_err).T