#This code requires three inputs: GC_name, mc_num and age
#This code outputs a list of chi2 for dm and reddening combination
import numpy as np
import pandas as pd
import os
from vorbin.voronoi_2d_binning import voronoi_2d_binning

class utiles:
	def check_file(self,path):
		if os.path.exists(path) == False:
			raise Exception("Cannot find inputfile at {}".format(path))
	
	def check_directories(self,path):
		if os.path.exists(path) == False:
			os.mkdir(path)
	
	#search for bin number for each data point
	def search_point_location_bc(self,x, y, xBar, yBar):
		bin_num = np.argmin(np.square(xBar - x) + np.square(yBar - y), axis = 1)
		return bin_num
		
	def writevorbin(self, xNode, yNode, bin_count_std, mc_num, age, path):
		#save the vorbin information
		bin_loc = np.vstack((xNode,yNode,bin_count_std)).T
		dp = pd.DataFrame(data=bin_loc, columns = ['xNode', 'yNode','bin_count_std'])
		path = "{}/bin_mc{}.age{}".format(path,mc_num,age)
		dp.to_csv(path, index=False)
	
	def search_vorbin(self, xNode, yNode, total_pt, vi, v):
		bin_count = np.zeros(len(xNode))
		Tb_size = self.Tb_size
		n_div = total_pt // Tb_size
		for i in range(n_div):
			bin_num = self.search_point_location_bc(v[i*Tb_size:(i+1)*Tb_size].reshape(Tb_size,1), vi[i*Tb_size:(i+1)*Tb_size].reshape(Tb_size,1), xNode, yNode)
			for j in range(Tb_size):
				bin_count[bin_num[j]] += 1
		#do the last bit
		len_last = total_pt - Tb_size * n_div
		if len_last != 0:
			bin_num = self.search_point_location_bc(v[-len_last:].reshape(len_last,1), vi[-len_last:].reshape(len_last,1), xNode, yNode)
			for j in range(len_last):
				bin_count[bin_num[j]] += 1
		#to avoid divde by 0
		for i in range(len(bin_count)):
			if bin_count[i] == 0:
				bin_count[i] += 1
		return bin_count
	
	#Divide the input data into three regions: MS, MSTO, and GB
	#required two cuts from the main
	def divide_data(self, df):
		df_MS = df[df['v'] > self.MSTO_cut]
		df_MSTO = df[(df['vi'] <= self.MSTO_cut) & (df['v'] >= self.GB_cut)]
		df_GB = df[df['vi'] < self.GB_cut]
		V_MS = df_MS['v'].values
		VI_MS = df_MS['vi'].values
		V_MSTO = df_MSTO['v'].values
		VI_MSTO = df_MSTO['vi'].values
		V_GB = df_GB['v'].values
		VI_GB = df_GB['vi'].values
		return V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB

	#generate vorbin for three different regions in the CMD: MS, MSTO, and GB
	#In default, we use 800 bins with 2:5:1 for three regions.
	def generate_vorbin(self,V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB,MS_bin_num=200,MSTO_bin_num=500,GB_bin_num=100,targetSN=10):
		#Do main sequence
		#Define number of stars used to generate vorbin based on TargetSN and bin_num
		fit_num = MS_bin_num*targetSN**2
		x = V_MS[:fit_num]
		y = VI_MS[:fit_num]
		signal = np.array([1]*fit_num)
		noise = signal
		#do vorbin
		_, x_gen_MS, y_gen_MS, _, _, _, _, _ = voronoi_2d_binning(x, y, signal, noise, targetSN, cvt=False, pixelsize=1, plot=False,quiet=True, sn_func=None, wvt=True)
		#Do main sequence turn off
		#Define number of stars used to generate vorbin based on TargetSN and bin_num
		fit_num = MSTO_bin_num*targetSN**2
		x = V_MSTO[:fit_num]
		y = VI_MSTO[:fit_num]
		signal = np.array([1]*fit_num)
		noise = signal
		#do vorbin
		_, x_gen_MSTO, y_gen_MSTO, _, _, _, _, _ = voronoi_2d_binning(x, y, signal, noise, targetSN, cvt=False, pixelsize=1, plot=False,quiet=True, sn_func=None, wvt=True)
		#Do giant branch
		#Define number of stars used to generate vorbin based on TargetSN and bin_num
		fit_num = GB_bin_num*targetSN**2
		x = V_GB[:fit_num]
		y = VI_GB[:fit_num]
		signal = np.array([1]*fit_num)
		noise = signal
		#do vorbin
		_, x_gen_GB, y_gen_GB, _, _, _, _, _ = voronoi_2d_binning(x, y, signal, noise, targetSN, cvt=False, pixelsize=1, plot=False,quiet=True, sn_func=None, wvt=True)
		#return location of bins
		x_gen = np.concatenate((x_gen_MS,x_gen_MSTO,x_gen_GB))
		y_gen = np.concatenate((y_gen_MS,y_gen_MSTO,y_gen_GB))
		return x_gen, y_gen



class chi2(utiles):

	def read_input(self,path):
		#read M92 observed data
		self.obs_data = pd.read_csv(path)
		#self.obs_size = len(self.obs_data)

	def main(self,write_vorbin,path,obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,dms,reds):
		#go through the search process
		width_coeff = (obs_v_max - obs_v_min)/(obs_vi_max - obs_vi_min)
		age = self.iso_age
		chi2 = []
		#read iso files
		dp = pd.read_csv("{}/mc{}.a{}".format(path,self.mc_num,age),sep='\s+',names=['vi','v'],skiprows=3)
		#filter out data points that is out of boundary
		df_cut = dp[(dp['vi'] < (obs_vi_max - red_max)) & (dp['vi'] > (obs_vi_min - red_min))& (dp['v'] < (obs_v_max - dm_max)) & (dp['v'] > (obs_v_min - dm_min))]
		total_pt = len(df_cut)
		#generate vorbin use the first 100000 data points
		V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB = self.divide_data(df_cut)
		x_gen, y_gen = self.generate_vorbin(V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB)
		#reduce memory usage for matrix operations
		XBar = np.float32(x_gen)
		YBar = np.float32(y_gen)
		v_32 = np.float32(df_cut['v'].values)
		vi_32 = np.float32(df_cut['vi'].values*width_coeff)
		#find standard bin count by search through all the theoretical data points
		bin_count_std = self.search_vorbin(XBar, YBar*width_coeff, total_pt, vi_32, v_32)
		#write vorbin infor if desired
		if write_vorbin == True:
			self.writevorbin(XBar, YBar, bin_count_std, self.mc_num, age,self.vorbin_path)
		#search through observed data
		for dm in dms:
			for red in reds:
				obs_cut = self.obs_data[(self.obs_data['vi'] - red < (obs_vi_max - red_max)) & (self.obs_data['vi'] - red > (obs_vi_min - red_min))& (self.obs_data['v'] - dm < (obs_v_max - dm_max)) & (self.obs_data['v'] - dm > (obs_v_min - dm_min))]
				obs_size = len(obs_cut)
				vi_32 = np.float32((obs_cut['vi'].values - red)*width_coeff)
				v_32 = np.float32(obs_cut['v'].values - dm)
				bin_count = self.search_vorbin(XBar, YBar*width_coeff, obs_size, vi_32, v_32)
				#calculate chi2
				chi2.append([age, dm, red, np.inner(np.divide(bin_count,bin_count_std/(total_pt/obs_size)) - 1, bin_count - bin_count_std/(total_pt/obs_size))])
		self.chi2 = chi2


	def writeout(self,path):
		#write chi2 to csv file
		dp = pd.DataFrame(data=self.chi2,columns=['age','dm','red','chi2'])
		path = "{}/chi2_a{}_mc{}".format(path,self.iso_age,self.mc_num)
		dp.to_csv(path)


	def __init__(self, GC_name, mc_num, iso_age, write_vorbin=False, Tb_size=30):
		#define boundaris
		obs_vi_max = 0.80
		obs_vi_min = 0.44
		obs_v_max = 19.28
		obs_v_min = 15.28
		#define distance modulus and reddening ranges
		dm_max = 13.90
		dm_min = 13.60
		red_max = 0.12
		red_min = 0.00
		dm_num = 31
		red_num = 7
		dms = np.linspace(dm_min,dm_max,dm_num)
		reds = np.linspace(red_min,red_max,red_num)
		#define other global variables
		self.mc_num = str(mc_num)
		self.iso_age = str(iso_age)
		self.Tb_size = Tb_size
		#define all the path for read and write
		obs_data_path = "/work2/08819/mying/{}/simulateCMD/{}_fitstars.dat".format(GC_name,GC_name)
		vorbin_path = "/work2/08819/mying/{}/vorbin".format(GC_name)
		self.vorbin_path = vorbin_path
		chi2_path = "/work2/08819/mying/{}/outchi2".format(GC_name)
		cmd_path = "/work2/08819/mying/{}/simulateCMD/outcmd".format(GC_name)
		#check those directories exist
		self.check_file(obs_data_path)
		self.check_directories(vorbin_path)
		self.check_directories(chi2_path)
		self.check_directories(cmd_path)
		#run code
		self.read_input(obs_data_path)
		self.main(write_vorbin,cmd_path,obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,dms,reds)        
		self.writeout(chi2_path)
		print("done mc{}".format(self.mc_num))

class resample(utiles):

	def read_input(self,obs_path,phot_path):
		#read M92 observed data
		self.obs_data = pd.read_csv(obs_path)
		#read AS test error
		self.dps_Ierr = []
		self.dps_Verr = []
		self.completeness_V = []
		self.completeness_I = []
		for i in range(80):
			self.dps_Ierr.append(pd.read_csv("{}\\Ierr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=3,names=['Ierr']))
			self.dps_Verr.append(pd.read_csv("{}\\Verr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=3,names=['Verr']))
			self.completeness_V.append(pd.read_csv("{}\\Verr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=1,names=["#","Npts","Radius","Mag","Completeness"],nrows=1)['Completeness'].values[0])
			self.completeness_I.append(pd.read_csv("{}\\Ierr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=1,names=["#","Npts","Radius","Mag","Completeness"],nrows=1)['Completeness'].values[0])

	def resample(self,k,path,write_cmd,obs_vi_max,obs_vi_min,obs_v_max,obs_v_min):
		sample_list = np.random.randint(0,len(self.obs_data),size=self.sample_pt)
		Ierr = np.zeros(self.sample_pt)
		Verr = np.zeros(self.sample_pt)
		Mask = [True]*self.sample_pt
		#print(self.obs_data['Ibin'] - 1)
		for i in range(self.sample_pt):
			Ierr[i] = self.dps_Ierr[self.obs_data['Ibin'].values[sample_list[i]]-1]['Ierr'].values[np.random.randint(0, high=len(self.dps_Ierr[self.obs_data['Ibin'].values[sample_list[i]]-1]))]
			Verr[i] = self.dps_Verr[self.obs_data['Vbin'].values[sample_list[i]]-1]['Verr'].values[np.random.randint(0, high=len(self.dps_Verr[self.obs_data['Vbin'].values[sample_list[i]]-1]))]
			#completeness test
			if np.random.rand() > (self.completeness_V[self.obs_data['Vbin'].values[sample_list[i]]-1]*self.completeness_I[self.obs_data['Ibin'].values[sample_list[i]]-1]):
				Mask[i] == False
		Vvega = self.obs_data['v'].values[sample_list] + Verr
		Ivega = self.obs_data['i'].values[sample_list] + Ierr
		VIvega = Vvega - Ivega
		Vvega = Vvega[Mask]
		VIvega= VIvega[Mask]
		data_resample = {'v':Vvega, 'vi':VIvega}
		dp = pd.DataFrame(data=data_resample)
		df_resample = dp[(dp['vi'] < (obs_vi_max)) & (dp['vi'] > (obs_vi_min))& (dp['v'] < (obs_v_max)) & (dp['v'] > (obs_v_min))]
		if write_cmd == True:
			path = "{}\\resample_{}".format(path,k)
			df_resample.to_csv(path,index=False)
		return df_resample

	def writeout(self,path,start,end):
		#write chi2 to csv file
		dp = pd.DataFrame(data=self.chi2,columns=['i','chi2'])
		path = "{}\\resample_chi2_{}_to_{}".format(path,start,end)
		dp.to_csv(path)


	def __init__(self, GC_name, start, end, MSTO_cut, GB_cut, write_vorbin=False, Tb_size=30,write_cmd=False, sample_pt=2000000):
		#define boundaris
		obs_vi_max = 0.791
		obs_vi_min = 0.473
		obs_v_max = 19.279
		obs_v_min = 15.297
		width_coeff = (obs_v_max - obs_v_min)/(obs_vi_max - obs_vi_min)
		#define all the path for read and write
		obs_data_path = "C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\GC_ages\\{}\\{}_fitstars_with_bins.dat".format(GC_name,GC_name)
		photometry_path = "C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\GC_ages\\{}\\inputfiles".format(GC_name)
		vorbin_path = "C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\GC_ages\\{}\\resample\\vorbin".format(GC_name)
		chi2_path = "C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\GC_ages\\{}\\resample\\outchi2".format(GC_name)
		cmd_path = "C:\\Users\\marti\\Desktop\\school work\\Dartmouth\\GC_ages\\{}\\resample\\cmd".format(GC_name)
		#check those directories exist
		self.check_file(obs_data_path)
		self.check_file(photometry_path)
		self.check_directories(vorbin_path)
		self.check_directories(chi2_path)
		self.check_directories(cmd_path)
		#assign other global variables
		self.Tb_size = Tb_size
		self.sample_pt = sample_pt
		#find both cuts from observational data
		self.MSTO_cut = MSTO_cut
		self.GB_cut = GB_cut
		#read obs data
		self.read_input(obs_data_path,photometry_path)
		#run resample
		self.chi2 = []
		for k in range(start, end):
			print("Starting {}th resample".format(k))
			df = self.resample(k,cmd_path,write_cmd,obs_vi_max,obs_vi_min,obs_v_max,obs_v_min)
			total_pt = len(df)
			V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB = self.divide_data(df)
			x_gen, y_gen = self.generate_vorbin(V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB)
			#reduce memory usage for matrix operations
			XBar = np.float32(x_gen)
			YBar = np.float32(y_gen)
			v_32 = np.float32(df['v'].values)
			vi_32 = np.float32(df['vi'].values*width_coeff)
			#find standard bin count by search through all the theoretical data points
			bin_count_std = self.search_vorbin(XBar, YBar*width_coeff, total_pt, vi_32, v_32)
			if write_vorbin == True:
				#no age for resample
				age = 0
				self.writevorbin(x_gen, y_gen, bin_count_std, k, age, vorbin_path)
			#fit observation data
			obs_size = len(self.obs_data)
			vi_32 = np.float32(self.obs_data['vi'].values*width_coeff)
			v_32 = np.float32(self.obs_data['v'].values)
			bin_count = self.search_vorbin(XBar, YBar*width_coeff, obs_size, vi_32, v_32)
			#calculate chi2
			self.chi2.append([k, np.inner(np.divide(bin_count,bin_count_std/(total_pt/obs_size)) - 1, bin_count - bin_count_std/(total_pt/obs_size))])  
		self.writeout(chi2_path,start,end)