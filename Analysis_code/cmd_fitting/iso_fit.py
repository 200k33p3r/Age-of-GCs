#simple isochrone fit code
import numpy as np
import pandas as pd
import os
from scipy.interpolate import CubicSpline
from skopt import gp_minimize
from global_var import define_range
from path_config import data_path
import sys

def read_iso(file_path):
    # Read the content of the file
    with open(file_path, 'r') as file:
        content = file.read()

    # Split the content into sub-files using two blank lines as the separator
    sub_files = content.split('\n\n\n')

    # Initialize a list to hold the dataframes with AGE(MYR) information
    iso_css = []
    ages = []

    for sub_file in sub_files:
        lines = sub_file.strip().split('\n')

        if len(lines) >= 4:  # Ensure there are enough lines to process
            local_header = ['EEP','Mass','v','vi']
            local_values = [line.split() for line in lines[3:]]  # Get the values from the rest lines

            # Create a dataframe for local information
            df = pd.DataFrame(local_values, columns=local_header).astype(float)

            # Extract the AGE(MYR) value from the second line
            age_value = lines[1].split()[3]
            ages.append(int(float(age_value)))

            # Append the dataframe to the list
            df_clean = df.drop_duplicates(subset=['v']).sort_values(by='v').copy()
            iso_css.append(CubicSpline(df_clean['v'].values, df_clean['vi'].values))
    return ages, iso_css

def dm_red_search(theta):
    dm,red = theta
    v_iso = np.linspace(V_SGB - dm - 2, V_SGB - dm + 2,1000)
    vi_iso = iso_cs(v_iso) + red
    return np.sum(np.square(vi_iso - vi_obs))

def simple_fit(mc_num, iso_path):
    file_path = iso_path + '/feh' + str(feh) + 'cmd.' + str(mc_num)
    try:
        ages, iso_css = read_iso(file_path)
        retval = []
        for i, age in enumerate(ages):
            global iso_cs
            iso_cs = iso_css[i]
            res = gp_minimize(dm_red_search,                  # the function to minimize
                      [(dm_min, dm_max), (red_min, red_max)],      # the bounds on each dimension of x
                      n_calls=500,         # the number of evaluations of f
                      n_random_starts=500,
                      verbose=False)
            chi2_fit, dm_fit, red_fit = res['fun'], res['x'][0], res['x'][1]
            retval.append([ age, dm_fit, red_fit, chi2_fit,mc_num])
        return pd.DataFrame(columns = ['age','dm','red','chi2','mc_num'],data =retval)
    except:
        pass

global GC_name, V_SGB, feh, dm_max, dm_min, red_max, red_min, vi_obs
GC_name = str(sys.argv[1])
mc_num = str(sys.argv[2])
V_SGB = 18.85
feh, dm_max, dm_min, red_max, red_min = define_range(GC_name)
iso_path = data_path + "{}/outiso".format(GC_name)
chi2_path = data_path + "{}/outchi2_iso".format(GC_name)
#read obs data
obs = pd.read_csv("{}_obs_iso.dat".format(GC_name))
vi_obs = obs['vi'].values
#calculdate chi2 and write
simple_fit(mc_num, iso_path).to_csv("{}/chi2_mc{}.csv".format(chi2_path, mc_num), index=False)