#This code is used to calculate the chi2 value
#using theoretical isochrones with radius, mass,
#and luminosity values to fit resampled eclipsing binaries.
import numpy as np
import pandas as pd
from rml_support import read_iso, calculate_chi2, check_file, writeout, resample_stars, read_candidates, read_eeps

def resample_run(iso_path, star_path, wrt_path, resample_num, resample_age, EEP_path = None, Consider_Lumi=True, Uniform_Resample=True):
    check_file(iso_path)
    check_file(star_path)
    Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err = read_candidates(star_path)
    mc_num, df_iso, _ = read_iso(iso_path)
    df_age = df_iso.loc[resample_age]
    if Uniform_Resample == False:
        EEPs_df = read_eeps(EEP_path)
        EEPs = EEPs_df.loc[(EEPs_df['mc_num'] == int(mc_num)) & (EEPs_df['Age'] == int(resample_age)), (EEPs_df.columns != 'mc_num') & (EEPs_df.columns != 'Age')].values
    else:
        EEPs = None
    Mass_star, Lumi_star, Rad_star = resample_stars(df_age, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, resample_num, EEPs = EEPs,Uniform_Resample= Uniform_Resample)
    chi2_data = calculate_chi2(df_age, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, Consider_Lumi = Consider_Lumi)
    dp = pd.DataFrame(data=chi2_data,columns=Names)
    dp['chi2_sum'] = np.sum(chi2_data,axis=1)
    writeout(dp, wrt_path)

def write_resample_star(iso_path, star_path, wrt_path, resample_age, EEP_path = None, Uniform_Resample=True):
    check_file(iso_path)
    check_file(star_path)
    Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err = read_candidates(star_path)
    mc_num, df_iso, _ = read_iso(iso_path)
    df_age = df_iso.loc[resample_age]
    if Uniform_Resample == False:
        EEPs_df = read_eeps(EEP_path)
        EEPs = EEPs_df.loc[(EEPs_df['mc_num'] == int(mc_num)) & (EEPs_df['Age'] == int(resample_age)), (EEPs_df.columns != 'mc_num') & (EEPs_df.columns != 'Age')].values
    else:
        EEPs = None
    Mass_star, Lumi_star, Rad_star = resample_stars(df_age, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, 1, EEPs = EEPs, Uniform_Resample= Uniform_Resample)
    dp = pd.read_csv(star_path)
    dp['M/Ms'] = Mass_star.flatten()
    dp['R/Rs'] = Rad_star.flatten()
    dp['L/Ls'] = Lumi_star.flatten()
    writeout(dp,wrt_path)

#def resample_multiple(iso_paths, star_path, wrt_path, )
