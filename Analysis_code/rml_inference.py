#This code is used to calculate the chi2 value
#using theoretical isochrones with radius, mass,
#and luminosity values to fit eclipsing binaries.

import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import os
from typing import Union
import numpy.typing as npt
FARRAY_1D = npt.NDArray[np.float64]

class utile:
    def check_file(self,path):
        if os.path.exists(path) == False:
            raise Exception("Cannot find inputfile at {}".format(path))
    
    def read_candidates(self,star_path):
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

    # def calculate_chi2(self, Mass,Lumi,Rad, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err):
    #     N_stars = len(Mass_star)
    #     N_eeps = len(Mass)
    #     chi2 = np.zeros((N_eeps,N_stars))
    #     #reshape for boardcasting
    #     Mass = Mass.reshape([1,N_eeps])
    #     Lumi = Lumi.reshape([1,N_eeps])
    #     Rad = Rad.reshape([1,N_eeps])
    #     #calculate difference
    #     Mass_diff = Mass - Mass_star
    #     Lumi_diff = Lumi - Lumi_star
    #     Rad_diff = Rad - Rad_star
    #     #find whether the difference is greater than 0 or not
    #     Mass_TF = Mass_diff > 0
    #     Lumi_TF = Lumi_diff > 0
    #     Rad_TF = Rad_diff > 0
    #     #define chi2
    #     Mass_chi2 = np.zeros((N_stars, N_eeps))
    #     Rad_chi2 = np.zeros((N_stars, N_eeps))
    #     Lumi_chi2 = np.zeros((N_stars, N_eeps))
    #     #for each star, evaluate whether the difference is positive or not. Then divide by the corresponding uncertainty
    #     for i in range(N_stars):
    #         Mass_chi2[i] = np.divide(Mass_diff[i],np.where(Mass_TF[i],Mass_star_p_err[i],Mass_star_m_err[i]))
    #         Lumi_chi2[i] = np.divide(Lumi_diff[i],np.where(Lumi_TF[i],Lumi_star_p_err[i],Lumi_star_m_err[i]))
    #         Rad_chi2[i] = np.divide(Rad_diff[i],np.where(Rad_TF[i],Rad_star_p_err[i],Rad_star_m_err[i]))
    #     #square the result
    #     Mass_chi2 = np.square(Mass_chi2)
    #     Lumi_chi2 = np.square(Lumi_chi2)
    #     Rad_chi2 = np.square(Rad_chi2)
    #     chi2 = Mass_chi2 + Lumi_chi2 + Rad_chi2
    #     indi_fit = np.min(chi2,axis=1)
    #     minchi2 = np.sum(indi_fit)
    #     return minchi2, indi_fit.tolist()
    
    def calculate_chi2(self, df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err):
        N_stars = len(Mass_star)
        #define the matrix value
        age_num = len(df_iso.index.get_level_values(0).unique())
        chi2_data = np.zeros((age_num, N_stars))
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
        return chi2_data

    def writeout(self,dp,wrt_path):
        dp.to_csv(wrt_path,index=False,mode='a',header=not os.path.exists(wrt_path))

    # def read_iso(self, iso_path, age='All',num_age=81):
    #     #read isochrone files with radius, mass and luminosity information
    #     self.check_file(iso_path)
    #     mc_num = int(iso_path[-5:])
    #     iso = open(iso_path, 'r')
    #     len_file = len(iso.readlines())
    #     iso.seek(0)
    #     len_idx = 1
    #     if len_file < 10:
    #         raise Exception("Empty file")
    #     else:
    #         if age == 'All':
    #             if num_age == 'Unknown':
    #                 num_age = 0
    #                 while len_idx <= len_file:
    #                     #look for how many ages it contains
    #                     iso.readline()
    #                     NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
    #                     num_age += 1
    #                     iso.readline()
    #                     len_idx += 3
    #                     for i in range(int(NPTS[1:])):
    #                         iso.readline()
    #                     len_idx += int(NPTS[1:])
    #                     iso.readline()
    #                     iso.readline()
    #                     len_idx += 2
    #                 #return to default
    #                 iso.seek(0)
    #                 len_idx = 1
    #             age_list = np.zeros(num_age)
    #             age_idx = 0
    #             while len_idx <= len_file:
    #                 #skip header
    #                 iso.readline()
    #                 NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
    #                 age_list[age_idx] = int(AGE[:-1])
    #                 age_idx += 1
    #                 Mass = np.zeros(int(NPTS[1:]))
    #                 Lumi = np.zeros(int(NPTS[1:]))
    #                 Rad = np.zeros(int(NPTS[1:]))
    #                 #skip header
    #                 iso.readline()
    #                 len_idx += 3
    #                 for i in range(int(NPTS[1:])):
    #                     EEP,MMs,LogG,LogTeff,LogLLs,LogRRs = iso.readline().split()
    #                     Mass[i] = float(MMs)
    #                     Lumi[i] = 10**(float(LogLLs))
    #                     Rad[i] = 10**(float(LogRRs))
    #                 len_idx += int(NPTS[1:])
    #                 iso.readline()
    #                 iso.readline()
    #                 len_idx += 2
    
    def read_iso(self, iso_path, max_eep=1000 ,age_list: Union[bool, FARRAY_1D] = np.linspace(8000,16000,81).astype(int),param_list=np.array(['Mass', 'Lumi', 'Rad'])):
        self.check_file(iso_path)
        mc_num = int(iso_path[-5:])
        iso = open(iso_path, 'r')
        len_file = len(iso.readlines())
        iso.seek(0)
        len_idx = 1
        if len_file < 10:
            raise Exception("Empty file")
        else:
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
            df = pd.DataFrame(data,index=arrays)
        return mc_num, df, age_list
    
    def resample_stars(self, df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err):
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


class run(utile):
    def main(self,iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err):
        mc_num, df_iso, age_list = self.read_iso(iso_path)
        chi2_data = self.calculate_chi2(df_iso, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err)
        dp = pd.DataFrame(data=chi2_data,columns=Names)
        dp['mc_num'] = np.full(len(dp),mc_num)
        dp['chi2'] = np.sum(chi2_data,axis=1)
        dp['Age'] = age_list
        return dp
    
    def __init__(self, iso_path, star_path,wrt_path):
        self.check_file(iso_path)
        self.check_file(star_path)
        Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err = self.read_candidates(star_path)
        dp = self.main(iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err)
        self.writeout(dp,wrt_path)

# class simulate(utile):
# #Add Monte Carlo method to simulate the chi2 values from theoretical isochrones.

#     def main(self, iso_path, candidates, num_shifts, mc_num):
#         num_stars = len(candidates)
#         self.check_file(iso_path)
#         iso = open(iso_path, 'r')
#         len_file = len(iso.readlines())
#         iso.seek(0)
#         len_idx = 1
#         head = np.array(['mc_num','Age','chi2','num_shifts'])
#         data = []
#         if len_file < 10:
#             raise Exception("Empty file")
#         else:
#             while len_idx <= len_file:
#                 #skip header
#                 iso.readline()
#                 NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
#                 Mass = []
#                 Lumi = []
#                 Rad = []
#                 #skip header
#                 iso.readline()
#                 len_idx += 3
#                 for i in range(int(NPTS[1:])):
#                     EEP,MMs,LogG,LogTeff,LogLLs,LogRRs = iso.readline().split()
#                     Mass.append(float(MMs))
#                     Lumi.append(10**(float(LogLLs)))
#                     Rad.append(10**(float(LogRRs)))
#                 Mass = np.array(Mass)
#                 Lumi = np.array(Lumi)
#                 Rad = np.array(Rad)
#                 #Pick only the q1 to q3 values to avoid potential issue
#                 NPT = len(Mass)
#                 # assume a uniform distribution along the isochone 
#                 # Can be updated to use some relation like PDMF
#                 for i in range(num_shifts):
#                     pick_idx = np.random.choice(range(int(NPT/4),NPT - int(NPT/4)), size=num_stars, replace=False)
#                     Mass_alter = Mass[pick_idx]
#                     Lumi_alter = Lumi[pick_idx]
#                     Rad_alter = Rad[pick_idx]
#                     sample_error= pd.DataFrame(columns = ['+dM/Ms','-dM/Ms','+dL/Ls','-dL/Ls','+dR/Rs','-dR/Rs'], data = np.abs(np.random.normal(scale=candidates[['+dM/Ms','-dM/Ms','+dL/Ls','-dL/Ls','+dR/Rs','-dR/Rs']].values)))
#                     sample_error[['-dM/Ms','-dL/Ls','-dR/Rs']] = -sample_error[['-dM/Ms','-dL/Ls','-dR/Rs']]
#                     candidates['M/Ms'] = Mass_alter + np.array([np.random.choice(sample_error[['+dM/Ms','-dM/Ms']].values[i]) for i in range(num_stars)])
#                     candidates['R/Rs'] = Rad_alter + np.array([np.random.choice(sample_error[['+dR/Rs','-dR/Rs']].values[i]) for i in range(num_stars)])
#                     candidates['L/Ls'] = Lumi_alter + np.array([np.random.choice(sample_error[['+dL/Ls','-dL/Ls']].values[i]) for i in range(num_stars)])
#                     minchi2, indi_fit = self.calculate_chi2(Mass,Lumi,Rad,candidates)
#                     data.append([mc_num,int(AGE[:-1]), minchi2, i])
#                 len_idx += int(NPTS[1:])
#                 iso.readline()
#                 iso.readline()
#                 len_idx += 2
#         dp = pd.DataFrame(data=data,columns=head)
#         return dp



#     def __init__(self, iso_path, star_path, wrt_path, num_shifts):
#         self.check_file(iso_path)
#         self.check_file(star_path)
#         candidates = self.read_candidates(star_path)
#         #assume the uncertainty from the candidates are "correct" and use them directly.
#         mc_num = int(iso_path[-5:])
#         dp = self.main(iso_path,candidates,num_shifts,mc_num)
#         self.writeout(dp,wrt_path)


# class simulate_iso(utile):
#     # extract stars from isochrones, add noise from candidates, rerun the analysis
    
#     def __init__(self, iso_path, star_path,wrt_path, resample_num):
#         self.check_file(iso_path)
#         self.check_file(star_path)
#         Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err = self.read_candidates(star_path)
#         mc_num, df_iso, age_list = self.read_iso(iso_path)
#         age_num = len(age_list)
#         # resample all the age in an isochrone file
#         retval = np.zeros((resample_num*age_num,3))
#         idx = 0
#         for resample_age in age_list:
#             for i in range(resample_num):
#                 df_iso_age = df_iso.loc[resample_age:resample_age+1].copy()
#                 Mass_star, Lumi_star, Rad_star = self.resample_stars(df_iso.loc[resample_age], Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err)
#                 chi2_data = self.calculate_chi2(df_iso_age, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err)
#                 chi2 = np.sum(chi2_data,axis=1)[0]
#                 retval[idx,0] = i
#                 retval[idx,1] = chi2
#                 retval[idx,2] = resample_age
#                 idx += 1
#         dp = pd.DataFrame(data=retval,columns=['num_sim', 'chi2', 'age'])
#         self.writeout(dp,wrt_path)