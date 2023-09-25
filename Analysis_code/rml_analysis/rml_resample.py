#This code is used to calculate the chi2 value
#using theoretical isochrones with radius, mass,
#and luminosity values to fit resampled eclipsing binaries.
import numpy as np
import pandas as pd
from rml_support import read_iso, calculate_chi2, check_file, writeout, resample_stars, read_candidates

def resample_run(iso_path, star_path, wrt_path, resample_num, resample_age, Consider_Lumi=True):
    check_file(iso_path)
    check_file(star_path)
    Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err = read_candidates(star_path)
    mc_num, df_iso, _ = read_iso(iso_path)
    df_age = df_iso.loc[resample_age]
    Mass_star, Lumi_star, Rad_star = resample_stars(df_age, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err,EEPs, resample_num)
    chi2_data = calculate_chi2(df_age, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, Consider_Lumi = Consider_Lumi)
    dp = pd.DataFrame(data=chi2_data,columns=Names)
    dp['chi2_sum'] = np.sum(chi2_data,axis=1)
    writeout(dp, wrt_path)

def resample_star(iso_path, star_path, wrt_path, resample_age):
    check_file(iso_path)
    check_file(star_path)
    Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err = read_candidates(star_path)
    mc_num, df_iso, _ = read_iso(iso_path)
    df_age = df_iso.loc[resample_age]
    Mass_star, Lumi_star, Rad_star = resample_stars(df_age, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err,EEPs, 1)
    dp = pd.read_csv(star_path)
    dp['M/Ms'] = Mass_star.flatten()
    dp['R/Rs'] = Rad_star.flatten()
    dp['L/Ls'] = Lumi_star.flatten()
    writeout(dp,wrt_path)
