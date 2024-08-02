from scipy.stats import qmc
import numpy as np
from scipy.stats import norm
import sys
import gs98getz
import pandas as pd

class genetatemcvar:

    def generate_sobol(self,m):
        #generate Sobol sequence to the correct dimension of space
        #num of vars in list + 2 (feh & afe)
        d = len(self.MC_var) + 2
        #add 1 to m to account for hard limit
        m += 1 
        sampler = qmc.Sobol(d)
        return sampler.random_base2(m)
    
    def write_out(self,outpath, out_df,feh):
        if feh < 0.0:
            filehead = '/varfeh-' + "{:03d}".format(int(np.abs(feh*100))) + '.'
        else:
            filehead = '/varfeh+' + "{:03d}".format(int(np.abs(feh*100))) + '.'
        for i, mc_num in enumerate(out_df['mc_num'].values):
            # Define the file to write to
            outfile = outpath + filehead + str(mc_num)
            # Open the file and write the data
            with open(outfile, 'w') as f:
                f.write("{:9.6f}     {:30s}\n".format(out_df['FeH'].values[i], 'FeH'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['He Abundance'].values[i], 'He Abundance'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['alphafe'].values[i], 'alphafe'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['CMIXLA'].values[i], 'CMIXLA'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['FGRY'].values[i], 'FGRY'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['FGRZ'].values[i], 'FGRZ'))
                f.write("{:2d}             {:30s}\n".format(out_df['TTAU'].values[i].astype(int), 'KTTAU'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['ALPHAE'].values[i], 'ALPHAE'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['ALPHAC'].values[i], 'ALPHAC'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['sstandard_1'].values[i], 'SSTANDARD(1) -- PP'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['sstandard_2'].values[i], 'SSTANDARD(2) -- He3+He3'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['sstandard_3'].values[i], 'SSTANDARD(3) -- He3+He4'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['sstandard_4'].values[i], 'SSTANDARD(4) -- P+C12'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['sstandard_5'].values[i], 'SSTANDARD(5) -- P+C13'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['sstandard_6'].values[i], 'SSTANDARD(6) -- P+N14'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['sstandard_7'].values[i], 'SSTANDARD(7) -- P+O16'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['alexcoef'].values[i], 'alexcoef'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['opalcoef2'].values[i], 'opalcoef2'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['talphacoef'].values[i], 'talphacoef'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['plascoef'].values[i], 'plascoef'))
                f.write("{:9.6f}     {:30s}\n".format(out_df['cocoef'].values[i], 'cocoef'))

    def __init__(self, outpath, mc_start, mc_end, feh = 0.0, sigma_feh = 0.0, afe = 0.0, sigma_afe = 0.0, min_bound = 0.05):
        #read MC variable lists
        MC_var = pd.read_csv('MC_var_ranges.csv')
        self.MC_var = MC_var.copy()
        #end exclusive
        Param_set_num = mc_end - mc_start
        m = np.floor(np.log2(Param_set_num))
        sobol_seq = self.generate_sobol(m)
        #check for hard limit
        bound_list = MC_var[MC_var['bound'] == True].index.tolist()
        #calculate the lower bound
        lower_bound = norm.cdf(min_bound,loc = MC_var.iloc[bound_list]['var2'].values,scale=MC_var.iloc[bound_list]['var1'].values)
        #compare values and create final mask
        masks = []
        for i in range(len(bound_list)):
            masks.append(sobol_seq[:,bound_list[i]] > lower_bound[i])
        final_mask = np.logical_and.reduce(masks)
        resample_runs = 0
        while len(sobol_seq[final_mask]) < Param_set_num:
            #resample
            if resample_runs > 5:
                #give up and increase the number of samples
                m += 1
            sobol_seq = self.generate_sobol(m)
            resample_runs += 1
            masks = []
            for i in range(len(bound_list)):
                masks.append(sobol_seq[:,bound_list[i]] > lower_bound[i])
            final_mask = np.logical_and.reduce(masks)
        final_sobol_seq = sobol_seq[final_mask][np.random.choice(np.arange(len(sobol_seq[final_mask])), size=Param_set_num, replace=False)]
        #now we convert the sobol sequence to actual output
        Out_MC_var = np.zeros(np.shape(final_sobol_seq))
        for i in range(len(MC_var)-2):
            if MC_var['dist'][i] == 'Uniform':
                Out_MC_var[:,i] = final_sobol_seq[:,i]*MC_var['var1'].iloc[i] + MC_var['var2'].iloc[i]
            elif MC_var['dist'][i] == 'Gaussian':
                Out_MC_var[:,i] = norm.ppf(final_sobol_seq[:,i])*MC_var['var1'].iloc[i]+MC_var['var2'].iloc[i]
            else:
                raise Exception('Unspecified distribution')
        #get feh
        Out_MC_var[:,-2] = norm.ppf(final_sobol_seq[:,-2])*sigma_feh + feh
        #get afe
        Out_MC_var[:,-1] = norm.ppf(final_sobol_seq[:,-2])*sigma_afe + afe
        #for [a/Fe], only have tables in steps of 0.2dex, so need to round to the nearest tabulated value
        for i in range(len(Out_MC_var[-1])):
            if Out_MC_var[-1,i] < -0.10:
                Out_MC_var[-1,i] = -0.20
            elif Out_MC_var[-1,i] < 0.10:
                Out_MC_var[-1,i] = 0.0
            elif Out_MC_var[-1,i] < 0.30:
                Out_MC_var[-1,i] = 0.20
            elif Out_MC_var[-1,i] < 0.50:
                Out_MC_var[-1,i] = 0.40
            elif Out_MC_var[-1,i] < 0.70:
                Out_MC_var[-1,i] = 0.60
            else:
                Out_MC_var[-1,i] = 0.80
        #we have 3 difference surface boundary conditions
        KATTU_idx = MC_var[MC_var['name']=='TTAU'].index.tolist()[0]
        for i in range(len(Out_MC_var[KATTU_idx])):
            if Out_MC_var[KATTU_idx,i] <= 1/3:
                Out_MC_var[KATTU_idx,i] = 0
            elif Out_MC_var[KATTU_idx,i] <= 2/3:
                Out_MC_var[KATTU_idx,i] = 1
            else:
                Out_MC_var[KATTU_idx,i] = 5
        #finally, get the correct Helium abundance
        z = np.zeros(len(Out_MC_var))
        for i in range(len(Out_MC_var)):
            sumz = gs98getz.gs98getz(Out_MC_var[i,-2], Out_MC_var[i,-1])
            z[i] = 0.244 + MC_var[MC_var['name']=='dyprim']['var1'].values*Out_MC_var[i,-3] + MC_var[MC_var['name']=='dyprim']['var2'].values + sumz*(MC_var[MC_var['name']=='dydz']['var1'].values*Out_MC_var[i,-4] + MC_var[MC_var['name']=='dydz']['var2'].values)
        Out_MC_var = np.hstack((Out_MC_var,z.reshape(len(z),1)))
        #write output as a pandas dataframe
        col_names = np.concatenate((MC_var['name'].values,np.array(['FeH','alphafe','He Abundance'])))
        out_df = pd.DataFrame(data=Out_MC_var, columns =col_names)
        out_df['mc_num'] = np.arange(mc_start,mc_end)
        #write outfile
        self.write_out(outpath, out_df, feh)