#This code is used to calculate the chi2 value
#using theoretical isochrones with radius, mass,
#and luminosity values to fit eclipsing binaries.

import numpy as np
import pandas as pd
import os

class utile:
    def check_file(self,path):
        if os.path.exists(path) == False:
            raise Exception("Cannot find inputfile at {}".format(path))
    
    def read_candidates(self,star_path):
        dp = pd.read_csv(star_path)
        N_stars = len(dp)
        Names = dp['Names'].values
        Mass_star = dp['M/Ms'].values.reshape(N_stars,1)
        Mass_star_p_err = dp['+dM/Ms'].values.reshape(N_stars,1)
        Mass_star_m_err = dp['-dM/Ms'].values.reshape(N_stars,1)
        Lumi_star = dp['L/Ls'].values.reshape(N_stars,1)
        Lumi_star_p_err = dp['+dL/Ls'].values.reshape(N_stars,1)
        Lumi_star_m_err = dp['-dL/Ls'].values.reshape(N_stars,1)
        Rad_star = dp['R/Rs'].values.reshape(N_stars,1)
        Rad_star_p_err = dp['+dR/Rs'].values.reshape(N_stars,1)
        Rad_star_m_err = dp['-dR/Rs'].values.reshape(N_stars,1)
        return Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err

    def calculate_chi2(self, Mass,Lumi,Rad, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err):
        N_stars = len(Mass_star)
        N_eeps = len(Mass)
        chi2 = np.zeros((N_eeps,N_stars))
        #reshape for boardcasting
        Mass = Mass.reshape([1,N_eeps])
        Lumi = Lumi.reshape([1,N_eeps])
        Rad = Rad.reshape([1,N_eeps])
        #calculate difference
        Mass_diff = Mass - Mass_star
        Lumi_diff = Lumi - Lumi_star
        Rad_diff = Rad - Rad_star
        #find whether the difference is greater than 0 or not
        Mass_TF = Mass_diff > 0
        Lumi_TF = Lumi_diff > 0
        Rad_TF = Rad_diff > 0
        #define chi2
        Mass_chi2 = np.zeros((N_stars, N_eeps))
        Rad_chi2 = np.zeros((N_stars, N_eeps))
        Lumi_chi2 = np.zeros((N_stars, N_eeps))
        #for each star, evaluate whether the difference is positive or not. Then divide by the corresponding uncertainty
        for i in range(N_stars):
            Mass_chi2[i] = np.divide(Mass_diff[i],np.where(Mass_TF[i],Mass_star_p_err[i],Mass_star_m_err[i]))
            Lumi_chi2[i] = np.divide(Lumi_diff[i],np.where(Lumi_TF[i],Lumi_star_p_err[i],Lumi_star_m_err[i]))
            Rad_chi2[i] = np.divide(Rad_diff[i],np.where(Rad_TF[i],Rad_star_p_err[i],Rad_star_m_err[i]))
        #square the result
        Mass_chi2 = np.square(Mass_chi2)
        Lumi_chi2 = np.square(Lumi_chi2)
        Rad_chi2 = np.square(Rad_chi2)
        chi2 = Mass_chi2 + Lumi_chi2 + Rad_chi2
        indi_fit = np.min(chi2,axis=1)
        minchi2 = np.sum(indi_fit)
        return minchi2, indi_fit.tolist()

    def writeout(self,dp,wrt_path):
        dp.to_csv(wrt_path,index=False,mode='a',header=not os.path.exists(wrt_path))


class run(utile):
    def main(self,iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err,mc_num):
        self.check_file(iso_path)
        iso = open(iso_path, 'r')
        len_file = len(iso.readlines())
        iso.seek(0)
        len_idx = 1
        head = np.concatenate((np.array(['mc_num','Age','chi2']),Names))
        data = []
        if len_file < 10:
            raise Exception("Empty file")
        else:
            while len_idx <= len_file:
                #skip header
                iso.readline()
                NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
                Mass = []
                Lumi = []
                Rad = []
                #skip header
                iso.readline()
                len_idx += 3
                for i in range(int(NPTS[1:])):
                    EEP,MMs,LogG,LogTeff,LogLLs,LogRRs = iso.readline().split()
                    Mass.append(float(MMs))
                    Lumi.append(10**(float(LogLLs)))
                    Rad.append(10**(float(LogRRs)))
                Mass = np.array(Mass)
                Lumi = np.array(Lumi)
                Rad = np.array(Rad)
                len_idx += int(NPTS[1:])
                iso.readline()
                iso.readline()
                len_idx += 2
                minchi2, indi_fit = self.calculate_chi2(Mass,Lumi,Rad,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err)
                data.append([mc_num,int(AGE[:-1]), minchi2] + indi_fit)
        dp = pd.DataFrame(data=data,columns=head)
        return dp
    
    def __init__(self, iso_path, star_path,wrt_path):
        self.check_file(iso_path)
        self.check_file(star_path)
        Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err = self.read_candidates(star_path)
        mc_num = int(iso_path[-5:])
        dp = self.main(iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, mc_num)
        self.writeout(dp,wrt_path)

class simulate(utile):
#Add Monte Carlo method to simulate the chi2 values from theoretical isochrones.

    def main(self, iso_path, candidates, num_shifts, mc_num):
        num_stars = len(candidates)
        self.check_file(iso_path)
        iso = open(iso_path, 'r')
        len_file = len(iso.readlines())
        iso.seek(0)
        len_idx = 1
        head = np.array(['mc_num','Age','chi2','num_shifts'])
        data = []
        if len_file < 10:
            raise Exception("Empty file")
        else:
            while len_idx <= len_file:
                #skip header
                iso.readline()
                NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
                Mass = []
                Lumi = []
                Rad = []
                #skip header
                iso.readline()
                len_idx += 3
                for i in range(int(NPTS[1:])):
                    EEP,MMs,LogG,LogTeff,LogLLs,LogRRs = iso.readline().split()
                    Mass.append(float(MMs))
                    Lumi.append(10**(float(LogLLs)))
                    Rad.append(10**(float(LogRRs)))
                Mass = np.array(Mass)
                Lumi = np.array(Lumi)
                Rad = np.array(Rad)
                #Pick only the q1 to q3 values to avoid potential issue
                NPT = len(Mass)
                # assume a uniform distribution along the isochone 
                # Can be updated to use some relation like PDMF
                for i in range(num_shifts):
                    pick_idx = np.random.choice(range(int(NPT/4),NPT - int(NPT/4)), size=num_stars, replace=False)
                    Mass_alter = Mass[pick_idx]
                    Lumi_alter = Lumi[pick_idx]
                    Rad_alter = Rad[pick_idx]
                    sample_error= pd.DataFrame(columns = ['+dM/Ms','-dM/Ms','+dL/Ls','-dL/Ls','+dR/Rs','-dR/Rs'], data = np.abs(np.random.normal(scale=candidates[['+dM/Ms','-dM/Ms','+dL/Ls','-dL/Ls','+dR/Rs','-dR/Rs']].values)))
                    sample_error[['-dM/Ms','-dL/Ls','-dR/Rs']] = -sample_error[['-dM/Ms','-dL/Ls','-dR/Rs']]
                    candidates['M/Ms'] = Mass_alter + np.array([np.random.choice(sample_error[['+dM/Ms','-dM/Ms']].values[i]) for i in range(num_stars)])
                    candidates['R/Rs'] = Rad_alter + np.array([np.random.choice(sample_error[['+dR/Rs','-dR/Rs']].values[i]) for i in range(num_stars)])
                    candidates['L/Ls'] = Lumi_alter + np.array([np.random.choice(sample_error[['+dL/Ls','-dL/Ls']].values[i]) for i in range(num_stars)])
                    minchi2, indi_fit = self.calculate_chi2(Mass,Lumi,Rad,candidates)
                    data.append([mc_num,int(AGE[:-1]), minchi2, i])
                len_idx += int(NPTS[1:])
                iso.readline()
                iso.readline()
                len_idx += 2
        dp = pd.DataFrame(data=data,columns=head)
        return dp



    def __init__(self, iso_path, star_path, wrt_path, num_shifts):
        self.check_file(iso_path)
        self.check_file(star_path)
        candidates = self.read_candidates(star_path)
        #assume the uncertainty from the candidates are "correct" and use them directly.
        mc_num = int(iso_path[-5:])
        dp = self.main(iso_path,candidates,num_shifts,mc_num)
        self.writeout(dp,wrt_path)


class simulate_iso(utile):
    # extract stars from isochrones, add noise from candidates, rerun the analysis
    def main(self, iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, mc_num, age, candidates):
        num_stars = len(Names)
        self.check_file(iso_path)
        iso = open(iso_path, 'r')
        len_file = len(iso.readlines())
        iso.seek(0)
        len_idx = 1
        if len_file < 10:
            raise Exception("Empty file")
        else:
            while len_idx <= len_file:
                #skip header
                iso.readline()
                NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
                Mass = []
                Lumi = []
                Rad = []
                #skip header
                iso.readline()
                len_idx += 3
                for i in range(int(NPTS[1:])):
                    EEP,MMs,LogG,LogTeff,LogLLs,LogRRs = iso.readline().split()
                    Mass.append(float(MMs))
                    Lumi.append(10**(float(LogLLs)))
                    Rad.append(10**(float(LogRRs)))
                Mass = np.array(Mass)
                Lumi = np.array(Lumi)
                Rad = np.array(Rad)
                #Pick only the q1 to q3 values to avoid potential issue
                NPT = len(Mass)
                # assume a uniform distribution along the isochone 
                # Can be updated to use some relation like PDMF
                if age == int(AGE[:-1]):
                    pick_idx = np.random.choice(range(int(NPT/4),NPT - int(NPT/4)), size=num_stars, replace=False)
                    Mass_alter = Mass[pick_idx]
                    Lumi_alter = Lumi[pick_idx]
                    Rad_alter = Rad[pick_idx]
                    sample_error= pd.DataFrame(columns = ['+dM/Ms','-dM/Ms','+dL/Ls','-dL/Ls','+dR/Rs','-dR/Rs'], data = np.abs(np.random.normal(scale=np.concatenate((Mass_star_p_err, Mass_star_m_err, Lumi_star_p_err, Lumi_star_m_err, Rad_star_p_err, Rad_star_m_err),axis=1))))
                    sample_error[['-dM/Ms','-dL/Ls','-dR/Rs']] = -sample_error[['-dM/Ms','-dL/Ls','-dR/Rs']]
                    candidates['M/Ms'] = Mass_alter + np.array([np.random.choice(sample_error[['+dM/Ms','-dM/Ms']].values[i]) for i in range(num_stars)])
                    candidates['R/Rs'] = Rad_alter + np.array([np.random.choice(sample_error[['+dR/Rs','-dR/Rs']].values[i]) for i in range(num_stars)])
                    candidates['L/Ls'] = Lumi_alter + np.array([np.random.choice(sample_error[['+dL/Ls','-dL/Ls']].values[i]) for i in range(num_stars)])
                    return candidates
                
                len_idx += int(NPTS[1:])
                iso.readline()
                iso.readline()
                len_idx += 2

    def __init__(self, iso_path, star_path, wrt_path, age):
        self.check_file(iso_path)
        self.check_file(star_path)
        candidates = pd.read_csv(star_path)
        Names, Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err = self.read_candidates(star_path)
        #assume the uncertainty from the candidates are "correct" and use them directly.
        mc_num = int(iso_path[-5:])
        dp = self.main(iso_path,Names,Mass_star, Mass_star_p_err, Mass_star_m_err, Lumi_star, Lumi_star_p_err, Lumi_star_m_err, Rad_star, Rad_star_p_err, Rad_star_m_err, mc_num, age, candidates)
        self.writeout(dp,wrt_path)