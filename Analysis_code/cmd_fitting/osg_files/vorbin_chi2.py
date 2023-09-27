#This code requires three inputs: GC_name, mc_num and age
#This code outputs a list of chi2 for dm and reddening combination
import numpy as np
import pandas as pd
import os
from vorbin.voronoi_2d_binning import voronoi_2d_binning
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.optimize import differential_evolution as DE

#Modify the path to be run on OSG

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
	def divide_data(self, df, read_track=False):
		if read_track==False:
			df_MS = df[df['v'] > self.MSTO_cut]
			df_MSTO = df[(df['vi'] <= self.MSTO_cut) & (df['v'] >= self.GB_cut)]
			df_GB = df[df['vi'] < self.GB_cut]
		else:
			MSTO_cut, GB_cut = read_track
			df_MS = df[df['v'] > MSTO_cut]
			df_MSTO = df[(df['v'] <= MSTO_cut) & (df['v'] >= GB_cut)]
			df_GB = df[df['v'] < GB_cut]
		V_MS = df_MS['v'].values
		VI_MS = df_MS['vi'].values*self.width_coeff
		V_MSTO = df_MSTO['v'].values
		VI_MSTO = df_MSTO['vi'].values*self.width_coeff
		V_GB = df_GB['v'].values
		VI_GB = df_GB['vi'].values*self.width_coeff
		return [V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB]
	
	#find MSTO and GB eeps when using isochrones to generate vorbin
	def find_two_eeps(self,path):
		iso = open("{}/feh{}.{}".format(path,self.feh,self.mc_num), 'r')
		len_file = len(iso.readlines())
		iso.seek(0)
		AGE_wanted = np.float(self.iso_age)
		if len_file < 10:
			raise Exception("Empty file")
		else:   
			#skip header
			iso.readline()
			NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
			while int(AGE[:-1]) != AGE_wanted:
				#skip header
				iso.readline()
				for i in range(int(NPTS[1:])):
					#skiplines
					iso.readline()
				#skiplines
				iso.readline()
				iso.readline()
				iso.readline()
				NPTS,MIXLEN,OVERSH,AGE,Y,Z,ZEFF,FeH,alphaFe = iso.readline().split()
			#skip header
			iso.readline()
			for i in range(int(NPTS[1:])):
				EEP,MASS,LogG,LogTeff,LogL,_,_,_,V,VI,F606W,F606WF814W = iso.readline().split()
				if int(EEP[:3]) == 474:
					MSTO_cut = float(F606W)
				elif int(EEP[:3]) == 530:
					GB_cut = float(F606W)
					return MSTO_cut, GB_cut

	#generate vorbin for three different regions in the CMD: MS, MSTO, and GB
	#In default, we use 800 bins with 2:5:1 for three regions.
	def generate_vorbin(self,Mags,MS_bin_num=200,MSTO_bin_num=500,GB_bin_num=100,UniSN=False, targetSN=10):
		if UniSN == True:
			V, VI = Mags
			#Define number of stars used to generate vorbin based on TargetSN and bin_num
			fit_num = 800*targetSN**2
			mask = np.random.choice(range(len(V)), size=fit_num, replace=False)
			x = V[mask]
			y = VI[mask]
			signal = np.array([1]*fit_num)
			noise = signal
			#do vorbin
			_, x_gen, y_gen, _, _, _, _, _ = voronoi_2d_binning(x, y, signal, noise, targetSN, cvt=False, pixelsize=1, plot=False,quiet=True, sn_func=None, wvt=True)
		else:
			V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB = Mags
			#Do main sequence
			#Define number of stars used to generate vorbin based on TargetSN and bin_num
			fit_num = MS_bin_num*targetSN**2
			mask = np.random.choice(range(len(V_MS)), size=fit_num, replace=False)
			x = V_MS[mask]
			y = VI_MS[mask]
			signal = np.array([1]*fit_num)
			noise = signal
			#do vorbin
			_, x_gen_MS, y_gen_MS, _, _, _, _, _ = voronoi_2d_binning(x, y, signal, noise, targetSN, cvt=False, pixelsize=1, plot=False,quiet=True, sn_func=None, wvt=True)
			#Do main sequence turn off
			#Define number of stars used to generate vorbin based on TargetSN and bin_num
			fit_num = MSTO_bin_num*targetSN**2
			mask = np.random.choice(range(len(V_MSTO)), size=fit_num, replace=False)
			x = V_MSTO[mask]
			y = VI_MSTO[mask]
			signal = np.array([1]*fit_num)
			noise = signal
			#do vorbin
			_, x_gen_MSTO, y_gen_MSTO, _, _, _, _, _ = voronoi_2d_binning(x, y, signal, noise, targetSN, cvt=False, pixelsize=1, plot=False,quiet=True, sn_func=None, wvt=True)
			#Do giant branch
			#Define number of stars used to generate vorbin based on TargetSN and bin_num
			fit_num = GB_bin_num*targetSN**2
			mask = np.random.choice(range(len(V_GB)), size=fit_num, replace=False)
			x = V_GB[mask]
			y = VI_GB[mask]
			signal = np.array([1]*fit_num)
			noise = signal
			#do vorbin
			_, x_gen_GB, y_gen_GB, _, _, _, _, _ = voronoi_2d_binning(x, y, signal, noise, targetSN, cvt=False, pixelsize=1, plot=False,quiet=True, sn_func=None, wvt=True)
			#return location of bins
			x_gen = np.concatenate((x_gen_MS,x_gen_MSTO,x_gen_GB))
			y_gen = np.concatenate((y_gen_MS,y_gen_MSTO,y_gen_GB))
		return x_gen, y_gen

	def find_fit_val(self,arr,x):
		for idx in range(len(arr)):
			if x < arr[idx]:
				return arr[idx], idx
		return arr[-1], idx

	#fine the empirical cdf from observational datapoints
	def empirical_cdf(self,v):
		v = np.sort(v)
		len_v = len(v)
		cdf = np.linspace(0,1,len_v)
		v,idx = np.unique(v, return_index=True)
		cdf = cdf[idx]
		return v, cdf
	
	def writeout_resample(self,path,start,end,chi2):
		#write chi2 to csv file
		dp = pd.DataFrame(data=chi2,columns=['i','chi2'])
		path = "{}/resample_chi2_{}_to_{}".format(path,start,end)
		dp.to_csv(path,index=False)
	
	#trim down the obseravtional data to the fiting region
	def obs_cut(self, dm, red):
		return self.obs_data[(self.obs_data[:,self.vi_idx] - red < (self.obs_vi_max)) & (self.obs_data[:,self.vi_idx] - red > (self.obs_vi_min))& (self.obs_data[:,self.v_idx] - dm < (self.obs_v_max)) & (self.obs_data[:,self.v_idx] - dm > (self.obs_v_min))]

	#find chi2/df
	# def dm_red_search(self, dm, red):
	#change for skopt
	def dm_red_search(self, theta):
		dm, red = theta
		new_obs = self.obs_cut(dm, red)
		obs_size = len(new_obs)
		bin_count = self.search_vorbin(self.XBar, self.YBar, obs_size, (new_obs[:,self.vi_idx] - red)*self.width_coeff, new_obs[:,self.v_idx] - dm)
		return (np.inner(np.divide(bin_count,self.bin_count_std/(self.total_pt/obs_size)) - 1, bin_count - self.bin_count_std/(self.total_pt/obs_size)))/(obs_size - 22)
	
class chi2(utiles):

	def read_input(self,path):
		#read M92 observed data
		obs_data = pd.read_csv(path)
		names = ['v','v_err','i','i_err','vi','vi_err','x','y']
		for i in range(len(names)):
			for j in range(len(obs_data.columns)):
				if names[i] == obs_data.columns[j]:
					setattr(self,names[i] + '_idx', j)
		self.obs_data = pd.read_csv(path).to_numpy()
		#self.obs_size = len(self.obs_data)

	def main(self,write_vorbin,path, dm_max, dm_min, red_max, red_min, iso_path,chi2_path,write_chi2_log=False,UniSN=False):
		age = self.iso_age
		#read cmd files
		cmd = pd.read_csv("{}/mc{}.a{}".format(path,self.mc_num,age),sep='\s+',names=['vi','v'],skiprows=3)
		#go through the search process
		self.obs_vi_max = max(cmd['vi'].values)
		self.obs_vi_min = min(cmd['vi'].values)
		self.obs_v_max = max(cmd['v'].values)
		self.obs_v_min = min(cmd['v'].values)
		self.width_coeff = (self.obs_v_max - self.obs_v_min)/(self.obs_vi_max - self.obs_vi_min)
		#check whether vorbin already existed or not
		if os.path.isfile("{}/bin_mc{}.age{}".format(self.vorbin_path,self.mc_num,age)) == True:
			vorbin = pd.read_csv("{}/bin_mc{}.age{}".format(self.vorbin_path,self.mc_num,age), skiprows=1,names=['x_gen', 'y_gen','bin_count_std'])
			self.XBar = np.float32(vorbin['x_gen'].values)
			self.YBar = np.float32(vorbin['y_gen'].values)
			self.bin_count_std = vorbin['bin_count_std'].values
			self.total_pt = np.sum(self.bin_count_std)
		else:
			if UniSN == False:
				#generate vorbin use the first 100000 data points
				V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB = self.divide_data(cmd,read_track=self.find_two_eeps(iso_path))
				x_gen, y_gen = self.generate_vorbin([V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB])
			else:
				#generate vorbin use the first 100000 data points
				x_gen, y_gen = self.generate_vorbin([cmd['v'].values,cmd['vi'].values*self.width_coeff], UniSN=True)
			#reduce memory usage for matrix operations
			self.XBar = np.float32(x_gen)
			self.YBar = np.float32(y_gen)
			v_32 = np.float32(cmd['v'].values)
			vi_32 = np.float32(cmd['vi'].values*self.width_coeff)
			#find standard bin count by search through all the theoretical data points
			self.total_pt = len(cmd)
			self.bin_count_std = self.search_vorbin(self.XBar, self.YBar, self.total_pt, vi_32, v_32)
			#write vorbin infor if desired
			if write_vorbin == True:
				self.writevorbin(self.XBar, self.YBar, self.bin_count_std, self.mc_num, age,self.vorbin_path)

		#implement differential Evolution algorithm
		bounds = [(dm_min, dm_max) ,(red_min,red_max)]
		res = DE(self.dm_red_search, bounds,tol=0.001,popsize=50)
		chi2_fit, dm_fit, red_fit = res['fun'], res['x'][0], res['x'][1]

		#define new bounds based on +/- 0.01 mag from the fit dm an +/- 0.01 mag in reddening
		dm_bound = np.linspace(dm_fit - 0.01, dm_fit + 0.01, 20)
		red_bound = np.linspace(red_fit - 0.01, red_fit + 0.01, 20)
		#run the grid search
		for dm in dm_bound:
			for red in red_bound:
				chi2 = self.dm_red_search([dm_fit,red_fit])
				if chi2 < chi2_fit:
					chi2_fit = chi2
					dm_fit = dm
					red_fit = red

		#write out the result
		retval = np.array([self.iso_age, dm_fit, red_fit, chi2_fit, len(self.obs_cut(dm_fit, red_fit))])
		#print(chi2)
		pd.DataFrame(retval).to_csv("{}/chi2_a{}_mc{}".format(chi2_path,self.iso_age,self.mc_num),header=None, index=None)

	def __init__(self, GC_name, mc_num, iso_age, UniSN=False, write_vorbin=False, Tb_size=30):
		#define distance modulus and reddening ranges
		if GC_name == 'M55':
			dm_max = 14.1
			dm_min = 13.8
			red_max = 0.15
			red_min = 0.08
		#define other global variables
		self.mc_num = str(mc_num)
		self.iso_age = str(iso_age)
		self.Tb_size = Tb_size
		if GC_name == 'M92':
			self.feh = 230
		elif GC_name == 'M55':
			self.feh = 190
		#define all the path for read and write
		cwd = os.getcwd()
		obs_data_path = cwd + "/{}_fitstars.dat".format(GC_name)
		vorbin_path = cwd + "/vorbin"
		self.vorbin_path = vorbin_path
		chi2_path = cwd + "/outchi2"
		cmd_path = cwd + "/outcmd"
		iso_path = cwd + "/outiso"
		#check those directories exist
		self.check_file(obs_data_path)
		self.check_directories(vorbin_path)
		self.check_directories(chi2_path)
		self.check_directories(cmd_path)
		self.check_directories(iso_path)
		#run code
		self.read_input(obs_data_path)
		self.main(write_vorbin,cmd_path, dm_max, dm_min, red_max, red_min,iso_path,chi2_path,UniSN=UniSN)
		print("done mc{}".format(self.mc_num))