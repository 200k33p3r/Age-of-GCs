#This file includes all the support functions for RML analysis

import warnings
import numpy as np
from numpy import newaxis
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
    return dp

def determine_EEP(df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, Consider_Lumi=True):
    N_stars = len(Mass_star)
    #define the matrix value
    age_list = df_iso.index.get_level_values(0).unique()
    age_num = len(age_list)
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
            indi_fit = np.argmin(chi2,axis=1)
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
            indi_fit = np.argmin(chi2,axis=1)
            chi2_data[:,i] = indi_fit
    return age_list, chi2_data.astype(int)

# def calculate_chi2(df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, Consider_Lumi=True):
#     #Case 1: 1 set of isochrones fitted to 1 set of stars
#     #Case 2: Mutiple set of isochrones fitted to 1 set of stars
#     #Case 3: 1 set of isochrones fitted to mutiple sets of stars
#     #Case 4: Mutiple set of isochrones fitted to mutiple sets of stars
    
#     #Get iso values
#     Mass_iso = df_iso.loc[(slice(None),'Mass'),:].copy().values
#     if Consider_Lumi == True:
#         Lumi_iso = df_iso.loc[(slice(None),'Lumi'),:].copy().values
#     Rad_iso = df_iso.loc[(slice(None),'Rad'),:].copy().values
#     #Determine which case are we evaluating
#     N_stars, num_sample, _ = np.shape(Mass_star)
#     N_iso, _ = np.shape(Mass_iso)
#     #define the matrix value
#     chi2_data = np.zeros((max(N_iso, num_sample), N_stars))
#     #Case 4
#     if num_sample > 1 and N_iso > 1:
#         raise Exception("Cannot fit multiple isochrones to mutiple sets of stars simultaniously")
#     else:
#         #calculate difference
#         Mass = Mass_iso - Mass_star
#         Rad = Rad_iso - Rad_star
#         #find whether the difference is greater than 0 or not
#         Mass_TF = Mass > 0
#         Rad_TF = Rad > 0
#         #for each star, evaluate whether the difference is positive or not. Then divide by the corresponding uncertainty
#         Mass_chi2 = np.array([np.divide(Mass[i],np.where(Mass_TF[i],Mass_star_p_err[i],Mass_star_m_err[i])) for i in range(N_stars)])
#         Rad_chi2 = np.array([np.divide(Rad[i],np.where(Rad_TF[i],Rad_star_p_err[i],Rad_star_m_err[i])) for i in range(N_stars)])
#         #square the result
#         Mass_chi2 = np.square(Mass_chi2)
#         Rad_chi2 = np.square(Rad_chi2)
#         if Consider_Lumi == True:
#             #Calculate Lumi only if Consider_Lumi = True
#             Lumi = Lumi_iso - Lumi_star
#             Lumi_TF = Lumi > 0
#             Lumi_chi2 = np.array([np.divide(Lumi[i],np.where(Lumi_TF[i],Lumi_star_p_err[i],Lumi_star_m_err[i])) for i in range(N_stars)])
#             Lumi_chi2 = np.square(Lumi_chi2)
#             chi2 = Mass_chi2 + Lumi_chi2 + Rad_chi2
#         else:
#             chi2 = Mass_chi2 + Rad_chi2
#         chi2_data = np.min(chi2,axis=2)
#     #Case 3: np.shape(chi2_data) = (num_star, num_resample)
#     #Case 2: np.shape(chi2_data) = (num_star, num_iso)
#     #Case 1: np.shape(chi2_data) = (num_star, 1)
#     return chi2_data.T

def calculate_chi2(candidates_df, iso_df, chi2_columns,target_ages):
    isochrones = [iso_df[iso_df['Age'] == age] for age in target_ages]
    candidate_pos_err = candidates_df[['+d' + name for name in chi2_columns]].values
    candidate_neg_err = candidates_df[['-d' + name for name in chi2_columns]].values
    Nstars = len(candidates_df)
    choice = np.array([candidate_neg_err,candidate_pos_err])
    min_chi2_list = np.zeros((len(target_ages), Nstars))
    for i, isochrone in enumerate(isochrones):
        Diff = candidates_df[chi2_columns].values.astype(float) - isochrone[chi2_columns].values[:,newaxis].astype(float)
        condition = np.where(Diff > 0, 1, 0)
        Errs = np.choose(condition,choice)
        chi2_indiv = np.sum(np.square(np.divide(Diff,Errs)),axis=2)
        indiv_idx = np.argmin(chi2_indiv,axis=0)
        for j in range(Nstars):
            min_chi2_list[i,j] = chi2_indiv[indiv_idx[j],j]
    #return the min chi2 value for individual stars with respect to individual 
    return min_chi2_list
        


def writeout(dp,wrt_path):
    dp.to_csv(wrt_path,index=False,mode='a',header=not os.path.exists(wrt_path))

# def read_iso(iso_path, iso_headers, output_format, max_eep=1500 ,age_list: Union[bool, FARRAY_1D] = np.linspace(8000,16000,81).astype(int),param_list=np.array(['Mass', 'Lumi', 'Rad'])):
#     check_file(iso_path)
#     mc_num = int(iso_path[-5:])
#     iso = open(iso_path, 'r')
#     len_file = len(iso.readlines())
#     iso.seek(0)
#     len_idx = 1
#     if len_file < 10:
#         iso.close()
#         raise Exception("Empty file")
#     else:
#         if output_format == 'DSEP_iso':
#             if age_list == 'Unknown':
#                 age_list = []
#                 while len_idx <= len_file:
#                     #look for how many ages it contains
#                     iso.readline()
#                     RAW_NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
#                     npts = int(RAW_NPTS[1:])
#                     age_list.append(int(AGE[:-1]))
#                     iso.readline()
#                     len_idx += 3
#                     for i in range(npts):
#                         iso.readline()
#                     len_idx += npts
#                     iso.readline()
#                     iso.readline()
#                     len_idx += 2
#                 #return to default
#                 iso.seek(0)
#                 len_idx = 1
#                 age_list = np.array(age_list)
#             #generate Mutiindex dataframe
#             num_age = len(age_list)
#             num_param = len(param_list)
#             arrays = [np.repeat(age_list,num_param),np.repeat(param_list.reshape(1,num_param),num_age, axis=0).reshape(num_age*num_param)]
#             data = np.full((num_age*num_param, max_eep),np.inf)
#             data_idx = 0
#             while len_idx <= len_file:
#                 #skip header
#                 iso.readline()
#                 RAW_NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
#                 npts = int(RAW_NPTS[1:])
#                 # Mass = np.zeros(npts)
#                 # Lumi = np.zeros(npts)
#                 # Rad = np.zeros(npts)
#                 #skip header
#                 iso.readline()
#                 len_idx += 3
#                 for i in range(npts):
#                     EEP,MMs,LogG,LogTeff,LogLLs,LogRRs = iso.readline().split()
#                     # Mass[i] = float(MMs)
#                     # Lumi[i] = float(LogLLs)
#                     # Rad[i] = float(LogRRs)
#                     data[data_idx,i] = float(MMs)
#                     data[data_idx+1,i] = 10**float(LogLLs)
#                     data[data_idx+2,i] = 10**float(LogRRs)
#                 data_idx += 3
#                 len_idx += npts
#                 iso.readline()
#                 iso.readline()
#                 len_idx += 2
#             iso.close()
#         else:
#             iso.close()
#             #generate Mutiindex dataframe
#             num_age = len(age_list)
#             num_param = len(param_list)
#             arrays = [np.repeat(age_list,num_param),np.repeat(param_list.reshape(1,num_param),num_age, axis=0).reshape(num_age*num_param)]
#             data = np.full((num_age*num_param, max_eep),np.inf)
#             #read in isochrone file
#             df = pd.read_csv(iso_path)
#             df['Age'] = (df['Age']/10**6).astype(int)
#             for i, age in enumerate(age_list):
#                 data[3*i,df[df['Age'] == age]['EEP']] = df[df['Age'] == age]['MMs']
#                 data[3*i+1,df[df['Age'] == age]['EEP']] = 10**(df[df['Age'] == age]['LogLLs'])
#                 data[3*i+2,df[df['Age'] == age]['EEP']] = 10**(df[df['Age'] == age]['LogRRs'])
#         df = pd.DataFrame(data,index=arrays)
#     return mc_num, df, age_list
def read_iso(iso_path, output_format):
    check_file(iso_path)
    mc_num = int(iso_path[-5:])
    iso = open(iso_path, 'r')
    len_file = len(iso.readlines())
    iso.seek(0)
    iso.close()
    if len_file < 10:
        raise Exception("Empty file")
    else:
        #give up on old iso format
        # if output_format == 'DSEP_iso':
        # #original DSED isochrone format
        #     df = pd.read_csv(iso_path,sep='\s+',header=None)
        #     start_indices = df[df[0] == '#EEP'].index
        #     end_indices = df[df[0] == '#NPTS'].index
        #     retval = []
        #     for i in range(len(start_indices)):
        #         if i < len(start_indices) - 1:
        #             sin_age_iso = df.iloc[start_indices[i] + 1:end_indices[i+1]].dropna(axis=1)
        #         else:
        #             sin_age_iso = df.iloc[start_indices[i] + 1:].dropna(axis=1)
        #         eeps = sin_age_iso[0].values
        #         #remove *
        #         for j in range(len(eeps)):
        #             if eeps[j][-1] == '*':
        #                 eeps[j] = eeps[j][:-1]
        #         sin_age_iso[0] = eeps.astype(int)
        #         sin_age_iso.columns = iso_header
        #         sin_age_iso['Age'] = int(float(df.iloc[end_indices[i]+1][3]))
        #         retval.append(sin_age_iso)
        #     iso_df = pd.concat(retval)
        #     iso_df[iso_header].astype(float)
        if output_format == 'iso_D':
        #new format
            iso_header=['EEP','Age','M/Ms','L/Ls','R/Rs']
            Take_exp = ['L/Ls','R/Rs']
            iso_df = pd.read_csv(iso_path,names=iso_header,skiprows=1)
        else:
            raise Exception("{} format is not supported".format(output_format))
        if Take_exp != None:
            iso_df[Take_exp] = 10**(iso_df[Take_exp].values)
    return mc_num, iso_df



# def resample_stars(df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, resample_num,EEPs=None,Uniform_Resample=True):
#     #given the information about an isochrone, and the observed binaries, we redraw stars from the theoretical isochrone with know uncertainty.
#     #If Uniform_Resample == True, we draw stars uniformly from the center 50% of the isochrones
#     #Otherwise, we choose resample stars from +/- 3 eeps from the given eep point
#     num_stars = len(Mass_star)
#     iso_Mass = df_iso.loc['Mass'].values
#     iso_Lumi = df_iso.loc['Lumi'].values
#     iso_Rad = df_iso.loc['Rad'].values
#     if Uniform_Resample == True:
#         iso_TF = iso_Mass < np.inf
#         iso_Mass = iso_Mass[iso_TF]
#         iso_Lumi = iso_Lumi[iso_TF]
#         iso_Rad = iso_Rad[iso_TF]
#         pick_idx = np.random.choice(range(int(len(iso_Mass)/4),len(iso_Mass) - int(len(iso_Mass)/4)), size=(num_stars, resample_num), replace=True)
#     else:
#         pick_idx = np.array([np.random.choice(range(EEPs[i]-3,EEPs[i]+3), size=resample_num, replace=True) for i in range(num_stars)])
#     Mass_alter = iso_Mass[pick_idx]
#     Lumi_alter = iso_Lumi[pick_idx]
#     Rad_alter = iso_Rad[pick_idx]
#     M_p_re = np.abs(np.random.normal(scale=np.repeat(Mass_star_p_err.reshape(num_stars,1),resample_num,axis=1)))
#     M_m_re = - np.abs(np.random.normal(scale=np.repeat(Mass_star_m_err.reshape(num_stars,1),resample_num,axis=1)))
#     M_TF = np.random.choice(a=[1, 0], size=(num_stars, resample_num))
#     M_err = np.multiply(M_p_re,M_TF) + np.multiply(M_m_re, (M_TF - 1)*-1)
#     R_p_re = np.abs(np.random.normal(scale=np.repeat(Rad_star_p_err.reshape(num_stars,1),resample_num,axis=1)))
#     R_m_re = - np.abs(np.random.normal(scale=np.repeat(Rad_star_m_err.reshape(num_stars,1),resample_num,axis=1)))
#     R_TF = np.random.choice(a=[1, 0], size=(num_stars, resample_num))
#     R_err = np.multiply(R_p_re,R_TF) + np.multiply(R_m_re, (R_TF - 1)*-1)
#     L_p_re = np.abs(np.random.normal(scale=np.repeat(Lumi_star_p_err.reshape(num_stars,1),resample_num,axis=1)))
#     L_m_re = - np.abs(np.random.normal(scale=np.repeat(Lumi_star_m_err.reshape(num_stars,1),resample_num,axis=1)))
#     L_TF = np.random.choice(a=[1, 0], size=(num_stars, resample_num))
#     L_err = np.multiply(L_p_re,L_TF) + np.multiply(L_m_re, (L_TF - 1)*-1)
#     #return mass lumi and rad for resampled stars
#     return (Mass_alter + M_err).reshape((num_stars, resample_num, 1)), (Lumi_alter + L_err).reshape((num_stars, resample_num, 1)), (Rad_alter + R_err).reshape((num_stars, resample_num, 1))

def resample_stars(df_iso, candidates_df, resample_header):
    #randomly pick from the middle 80% of eeps in isochrones
    resample_idx = np.random.choice(range(int(0.1*len(df_iso)),int(0.9*len(df_iso))), size=len(candidates_df))
    #read in uncertainty
    candidate_pos_err = candidates_df[['+d' + name for name in resample_header]].values
    candidate_neg_err = candidates_df[['-d' + name for name in resample_header]].values
    #choose which size of curve to shift to
    err_idx = np.random.choice(a=[0,1],size=np.shape(candidate_pos_err))
    #take gaussian of uncertainty and add back to the selected eep point
    return candidates_df[resample_header].iloc[resample_idx].values + np.multiply(np.abs(np.random.normal(scale=candidate_pos_err)), err_idx) - np.multiply(np.abs(np.random.normal(scale=candidate_neg_err)), 1 - err_idx)

def read_eeps(eep_path):
    check_file(eep_path)
    df = pd.read_csv(eep_path,dtype=np.int32)
    return df