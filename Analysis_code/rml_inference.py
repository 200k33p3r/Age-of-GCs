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
        return dp

    def calculate_chi2(self, Mass,Lumi,Rad, candidates):
        N_stars = len(candidates)
        N_eeps = len(Mass)
        chi2 = np.zeros((N_eeps,N_stars))
        for i in range(N_stars):
            star_vals = candidates.iloc[i].values[1:].astype(float)
            for j in range(N_eeps):
                #fit mass
                if Mass[j] >= star_vals[0]:
                    chi2[j,i] += (Mass[j] - star_vals[0])**2/(star_vals[1]**2)
                else:
                    chi2[j,i] += (Mass[j] - star_vals[0])**2/(star_vals[2]**2)
                #fit luminosity
                if Lumi[j] >= star_vals[3]:
                    chi2[j,i] += (Lumi[j] - star_vals[3])**2/(star_vals[4]**2)
                else:
                    chi2[j,i] += (Lumi[j] - star_vals[3])**2/(star_vals[5]**2)
                #fit radius
                if Rad[j] >= star_vals[6]:
                    chi2[j,i] += (Rad[j] - star_vals[6])**2/(star_vals[7]**2)
                else:
                    chi2[j,i] += (Rad[j] - star_vals[6])**2/(star_vals[8]**2)
        indi_fit = np.min(chi2,axis=0)
        minchi2 = np.sum(indi_fit)
        return minchi2, indi_fit.tolist()

    def writeout(self,dp,wrt_path):
        dp.to_csv(wrt_path,index=False,mode='a',header=not os.path.exists(wrt_path))


class run(utile):
    def main(self,iso_path,candidates,mc_num):
        self.check_file(iso_path)
        iso = open(iso_path, 'r')
        len_file = len(iso.readlines())
        iso.seek(0)
        len_idx = 1
        head = np.concatenate((np.array(['mc_num','Age','chi2']),candidates['Names'].values))
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
                minchi2, indi_fit = self.calculate_chi2(Mass,Lumi,Rad,candidates)
                data.append([mc_num,int(AGE[:-1]), minchi2] + indi_fit)
        dp = pd.DataFrame(data=data,columns=head)
        return dp
    
    def __init__(self, iso_path, star_path,wrt_path):
        self.check_file(iso_path)
        self.check_file(star_path)
        candidates = self.read_candidates(star_path)
        mc_num = int(iso_path[-5:])
        dp = self.main(iso_path,candidates,mc_num)
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
    def main(self, iso_path, candidates, mc_num, age):
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
                if age == int(AGE[:-1]):
                    pick_idx = np.random.choice(range(int(NPT/4),NPT - int(NPT/4)), size=num_stars, replace=False)
                    Mass_alter = Mass[pick_idx]
                    Lumi_alter = Lumi[pick_idx]
                    Rad_alter = Rad[pick_idx]
                    sample_error= pd.DataFrame(columns = ['+dM/Ms','-dM/Ms','+dL/Ls','-dL/Ls','+dR/Rs','-dR/Rs'], data = np.abs(np.random.normal(scale=candidates[['+dM/Ms','-dM/Ms','+dL/Ls','-dL/Ls','+dR/Rs','-dR/Rs']].values)))
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
        candidates = self.read_candidates(star_path)
        #assume the uncertainty from the candidates are "correct" and use them directly.
        mc_num = int(iso_path[-5:])
        dp = self.main(iso_path,candidates,mc_num,age)
        self.writeout(dp,wrt_path)