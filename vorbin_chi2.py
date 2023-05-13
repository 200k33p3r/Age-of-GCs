#This code requires two inputs: mc_num and age
#This code outputs a list of chi2 for dm and reddening combination
import numpy as np
import pandas as pd
import sys
import os
from vorbin.voronoi_2d_binning import voronoi_2d_binning
import time

class chi2:

    def read_input(self):
        #read M92 observed data
        self.obs_data = pd.read_csv('/work2/08819/mying/M92/simulateCMD_ref/M92_fit_cut.dat')
        #self.obs_size = len(self.obs_data)

    #search for bin number for each data point
    def search_point_location_bc(self,x, y, xBar, yBar):
        bin_num = np.argmin(np.square(xBar - x) + np.square(yBar - y), axis = 1)
        return bin_num
        
    def writevorbin(self, xNode, yNode, bin_count_std, age):
        #save the vorbin information
        bin_loc = np.vstack((xNode,yNode,bin_count_std)).T
        dp = pd.DataFrame(data=bin_loc, columns = ['xNode', 'yNode','bin_count_std'])
        path = "/work2/08819/mying/M92/outchi2/bin_mc{}.age{}".format(self.mc_num,age)
        dp.to_csv(path, index=False)
    
    def search_vorbin(self, xNode, yNode, total_pt, vi, v):
        bin_count = np.zeros(len(xNode))
        Tb_size = self.Tb_size
        n_div = total_pt // Tb_size
        for i in range(n_div):
            bin_num = self.search_point_location_bc(vi[i*Tb_size:(i+1)*Tb_size].reshape(Tb_size,1), v[i*Tb_size:(i+1)*Tb_size].reshape(Tb_size,1), xNode, yNode)
            for j in range(Tb_size):
                bin_count[bin_num[j]] += 1
        #do the last bit
        len_last = total_pt - Tb_size * n_div
        if len_last != 0:
            bin_num = self.search_point_location_bc(vi[-len_last:].reshape(len_last,1), v[-len_last:].reshape(len_last,1), xNode, yNode)
            for j in range(len_last):
                bin_count[bin_num[j]] += 1
        #to avoid divde by 0
        for i in range(len(bin_count)):
            if bin_count[i] == 0:
                bin_count[i] += 1
        return bin_count

    def main(self):
        #go through the search process
        self.Tb_size = 30
        obs_vi_max = 0.721000016
        obs_vi_min = 0.351000011
        obs_v_max = 18.9250000
        obs_v_min = 16.9250000
        self.width_coeff = (obs_v_max - obs_v_min)/(obs_vi_max - obs_vi_min)
        dm_max = 14.82
        dm_min = 14.62
        red_max = 0.05
        red_min = 0.0
        vorbin_pt = 20000
        dms = np.linspace(dm_min,dm_max,21)
        reds = np.linspace(red_min,red_max,6)
        age = self.iso_age
        chi2 = []
        #read iso files
        dp = pd.read_csv("/work2/08819/mying/M92/simulateCMD/outcmd/mc{}.a{}".format(self.mc_num,age),sep='\s+',names=['vi','v'],skiprows=3)
        #filter out data points that is out of boundary
        dp = dp[(dp['vi'] < (obs_vi_max - red_max)) & (dp['vi'] > (obs_vi_min - red_min))& (dp['v'] < (obs_v_max - dm_max)) & (dp['v'] > (obs_v_min - dm_min))]
        total_pt = len(dp)
        #generate vorbin use the first 100000 data points
        y = np.float32(dp['v'].values[:vorbin_pt])
        x = np.float32(dp['vi'].values[:vorbin_pt] * self.width_coeff)
        signal = np.array([1]*vorbin_pt)
        noise = np.array([1]*vorbin_pt)
        targetSN = 5
        binNum, x_gen, y_gen, x_bar, y_bar, sn, nPixels, scale = voronoi_2d_binning(x, y, signal, noise, target_sn, cvt=False, pixelsize=1, plot=False,quiet=True, sn_func=None, wvt=True)
        #reduce memory usage for matrix operations
        XBar = np.float32(xBar)
        YBar = np.float32(yBar)
        #find standard bin count by search through all the theoretical data points
        bin_count_std = self.search_vorbin(XBar, YBar, len(dp), np.float32(dp['vi'].values*self.width_coeff), np.float32(dp['v'].values))
        self.writevorbin(XBar, YBar, bin_count_std, age)
        #search through observed data
        for dm in dms:
            for red in reds:
                bin_count = np.zeros(len(xBar))
                dp = self.obs_data[(self.obs_data['vi'] - red < (obs_vi_max - red_max)) & (self.obs_data['vi'] - red > (obs_vi_min - red_min))& (self.obs_data['v'] - dm < (obs_v_max - dm_max)) & (self.obs_data['v'] - dm > (obs_v_min - dm_min))]
                obs_size = len(dp)
                bin_count = self.search_vorbin(XBar, YBar, len(dp), np.float32((dp['vi'].values - red)*self.width_coeff), np.float32(dp['v'].values - dm))
                #calculate chi2
                chi2.append([age, dm, red, np.inner(np.divide(bin_count,bin_count_std/(total_pt/obs_size)) - 1, bin_count - bin_count_std/(total_pt/obs_size))])
        self.chi2 = chi2


    def writeout(self):
        #write chi2 to csv file
        dp = pd.DataFrame(data=self.chi2,columns=['age','dm','red','chi2'])
        path = "/work2/08819/mying/M92/outchi2/chi2_a{}_mc{}".format(self.iso_age,self.mc_num)
        dp.to_csv(path)


    def __init__(self, mc_num, iso_age):
        self.mc_num = str(mc_num)
        self.iso_age = str(iso_age)
        self.read_input()
        self.main()        
        self.writeout()
        print("done mc{}".format(self.mc_num))