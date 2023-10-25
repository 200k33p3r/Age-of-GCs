#This code is used to calculate the chi2 value
#using theoretical isochrones with radius, mass,
#and luminosity values to fit eclipsing binaries.

import numpy as np
import pandas as pd
from rml_support import read_iso, calculate_chi2, check_file, writeout, read_candidates, determine_EEP

# def main(iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, Consider_Lumi,iso_D):
#     mc_num, df_iso, age_list = read_iso(iso_path,iso_D=iso_D)
#     chi2_data = calculate_chi2(df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err,Consider_Lumi=Consider_Lumi)
#     dp = pd.DataFrame(data=chi2_data,columns=Names)
#     dp['mc_num'] = np.full(len(dp),mc_num)
#     dp['chi2'] = np.sum(chi2_data,axis=1)
#     dp['Age'] = age_list
#     return dp

# def run_analysis(iso_path, star_path,wrt_path,iso_format='iso_D', iso_header=['EEP','Age','M/Ms','L/Ls','R/Rs']):
#     name_lists = ['Names','M/Ms','+dM/Ms','-dM/Ms','L/Ls','+dL/Ls','-dL/Ls','R/Rs','+dR/Rs','-dR/Rs']
#     check_file(iso_path)
#     check_file(star_path)
#     dp = read_candidates(star_path,name_lists)
#     Names = dp['Names'].values
#     Mass_star = dp['M/Ms'].values
#     Mass_star_p_err = dp['+dM/Ms'].values
#     Mass_star_m_err = dp['-dM/Ms'].values
#     Lumi_star = dp['L/Ls'].values
#     Lumi_star_p_err = dp['+dL/Ls'].values
#     Lumi_star_m_err = dp['-dL/Ls'].values
#     Rad_star = dp['R/Rs'].values
#     Rad_star_p_err = dp['+dR/Rs'].values
#     Rad_star_m_err = dp['-dR/Rs'].values
#     dp = main(iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err)
#     writeout(dp,wrt_path)

def find_eeps(iso_paths, star_path, wrt_path, Consider_Lumi=True,iso_D=True):
    name_lists = ['Names','M/Ms','+dM/Ms','-dM/Ms','L/Ls','+dL/Ls','-dL/Ls','R/Rs','+dR/Rs','-dR/Rs']
    check_file(star_path)
    dp = read_candidates(star_path,name_lists)
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
    for iso_path in iso_paths:
        check_file(iso_path)
        mc_num, df_iso, age_list = read_iso(iso_path,iso_D=iso_D)
        age_list, chi2_data = determine_EEP(df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err,Consider_Lumi=Consider_Lumi)
        dp = pd.DataFrame(data=chi2_data,columns=Names)
        dp['Age'] = age_list
        dp['mc_num'] = [mc_num]*len(dp)
        writeout(dp,wrt_path)

def fit_rml(iso_path, star_path, wrt_path, iso_format='iso_D',metric_headers=['R/Rs','M/Ms','L/Ls'],target_ages=np.linspace(8000,16000,81).astype(int)):
    check_file(iso_path)
    check_file(star_path)
    #the default format is the new isochrone format
    if iso_format == 'iso_D':
        iso_header=['EEP','Age','M/Ms','L/Ls','R/Rs']
        Take_exp = ['L/Ls','R/Rs']
    elif iso_format == 'DSEP_iso':
        iso_header=['EEP','M/Ms','LogG','LogTeff','L/Ls','R/Rs']
        Take_exp = ['L/Ls','R/Rs']
    #read candidates
    Obs_stars = read_candidates(star_path)
    #read in isochrones
    mc_num, df_iso = read_iso(iso_path,iso_header,iso_format,Take_exp=Take_exp)
    #calculate chi2
    chi2_data = calculate_chi2(Obs_stars, df_iso, metric_headers,target_ages)
    dp = pd.DataFrame(data=chi2_data,columns=Obs_stars['Names'].values)
    dp['mc_num'] = np.full(len(dp),mc_num)
    dp['chi2'] = np.sum(chi2_data,axis=1)
    dp['Age'] = target_ages
    writeout(dp,wrt_path)

#calculate the chi2 value on the cmd plane for individual stars
def fit_cmd(iso_path, star_path, wrt_path, iso_format='DSEP_iso',metric_headers=['F606W','F606W-F814W'],target_ages=np.linspace(8000,16000,81).astype(int)):
    check_file(iso_path)
    check_file(star_path)
    #the default format is the old isochrone format
    if iso_format == 'iso_D':
        iso_header=['EEP','M/Ms','F606W','F606W-F814W']
        Take_exp=None
    #read candidates
    Obs_stars = read_candidates(star_path)
    #read in isochrones
    mc_num, df_iso = read_iso(iso_path,iso_header,iso_format,Take_exp=Take_exp)
    #calculate chi2
    chi2_data = calculate_chi2(Obs_stars, df_iso, metric_headers,target_ages)
    dp = pd.DataFrame(data=chi2_data,columns=Obs_stars['Names'].values)
    dp['mc_num'] = np.full(len(dp),mc_num)
    dp['chi2'] = np.sum(chi2_data,axis=1)
    dp['Age'] = target_ages
    writeout(dp,wrt_path)