#This code is used to calculate the chi2 value
#using theoretical isochrones with radius, mass,
#and luminosity values to fit eclipsing binaries.

import numpy as np
import pandas as pd
from rml_support import read_iso, calculate_chi2, check_file, writeout, read_candidates

def main(iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, Consider_Lumi,iso_D):
    mc_num, df_iso, age_list = read_iso(iso_path,iso_D=iso_D)
    chi2_data = calculate_chi2(df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err,Consider_Lumi=Consider_Lumi)
    dp = pd.DataFrame(data=chi2_data,columns=Names)
    dp['mc_num'] = np.full(len(dp),mc_num)
    dp['chi2'] = np.sum(chi2_data,axis=1)
    dp['Age'] = age_list
    return dp

def run_analysis(iso_path, star_path,wrt_path,Consider_Lumi=True,iso_D=True):
    check_file(iso_path)
    check_file(star_path)
    Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err = read_candidates(star_path)
    dp = main(iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err,Consider_Lumi,iso_D)
    writeout(dp,wrt_path)