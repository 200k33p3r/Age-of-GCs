#This code requires three inputs: GC_name, mc_num and age
#This code outputs a list of chi2 for dm and reddening combination
import numpy as np
import pandas as pd
import os
from vorbin.voronoi_2d_binning import voronoi_2d_binning
import subprocess
from scipy.interpolate import interp1d

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
			df_MSTO = df[(df['vi'] <= MSTO_cut) & (df['v'] >= GB_cut)]
			df_GB = df[df['vi'] < GB_cut]
		V_MS = df_MS['v'].values
		VI_MS = df_MS['vi'].values
		V_MSTO = df_MSTO['v'].values
		VI_MSTO = df_MSTO['vi'].values
		V_GB = df_GB['v'].values
		VI_GB = df_GB['vi'].values
		return V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB
	
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
	
	def writeout_resample(self,path,start,end):
		#write chi2 to csv file
		dp = pd.DataFrame(data=self.chi2,columns=['i','chi2'])
		path = "{}/resample_chi2_{}_to_{}".format(path,start,end)
		dp.to_csv(path)

class chi2(utiles):

	def read_input(self,path):
		#read M92 observed data
		self.obs_data = pd.read_csv(path)
		#self.obs_size = len(self.obs_data)

	def main(self,write_vorbin,path,obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,dms,reds,iso_path):
		#go through the search process
		width_coeff = (obs_v_max - obs_v_min)/(obs_vi_max - obs_vi_min)
		age = self.iso_age
		chi2 = []
		#read iso files
		dp = pd.read_csv("{}/mc{}.a{}".format(path,self.mc_num,age),sep='\s+',names=['vi','v'],skiprows=3)
		#filter out data points that is out of boundary
		df_cut = dp[(dp['vi'] < (obs_vi_max - reds[-1])) & (dp['vi'] > (obs_vi_min - reds[0]))& (dp['v'] < (obs_v_max - dms[-1])) & (dp['v'] > (obs_v_min - dms[0]))]
		total_pt = len(df_cut)
		#generate vorbin use the first 100000 data points

		V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB = self.divide_data(df_cut,read_track=self.find_two_eeps(iso_path))
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
				obs_cut = self.obs_data[(self.obs_data['vi'] - red < (obs_vi_max - reds[-1])) & (self.obs_data['vi'] - red > (obs_vi_min - reds[0]))& (self.obs_data['v'] - dm < (obs_v_max - dms[-1])) & (self.obs_data['v'] - dm > (obs_v_min - dms[0]))]
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
		if GC_name == 'M55':
			obs_vi_max = 0.792
			obs_vi_min = 0.462
			obs_v_max = 19.28
			obs_v_min = 15.296
		#define distance modulus and reddening ranges
		if GC_name == 'M55':
			dm_max = 14.00
			dm_min = 13.50
			red_max = 0.15
			red_min = 0.00
			dm_num = 51
			red_num = 16
		dms = np.linspace(dm_min,dm_max,dm_num)
		reds = np.linspace(red_min,red_max,red_num)
		#define other global variables
		self.mc_num = str(mc_num)
		self.iso_age = str(iso_age)
		self.Tb_size = Tb_size
		if GC_name == 'M92':
			self.feh = 230
		elif GC_name == 'M55':
			self.feh = 190
		#define all the path for read and write
		obs_data_path = "/work2/08819/mying/{}/simulateCMD/{}_fitstars.dat".format(GC_name,GC_name)
		vorbin_path = "/work2/08819/mying/{}/vorbin".format(GC_name)
		self.vorbin_path = vorbin_path
		chi2_path = "/work2/08819/mying/{}/outchi2".format(GC_name)
		cmd_path = "/work2/08819/mying/{}/simulateCMD/outcmd".format(GC_name)
		iso_path = "/work2/08819/mying/{}/outiso".format(GC_name)
		#check those directories exist
		self.check_file(obs_data_path)
		self.check_directories(vorbin_path)
		self.check_directories(chi2_path)
		self.check_directories(cmd_path)
		self.check_directories(iso_path)
		#run code
		self.read_input(obs_data_path)
		self.main(write_vorbin,cmd_path,obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,dms,reds,iso_path)        
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
		self.writeout_resample(chi2_path,start,end)

#this is similar to resample class but utilizing the fiducial isochrone generated from fidanka
class resample_fidanka(utiles):
	def read_input(self,obs_path,phot_path,fiducial_path):
		#read M92 observed data
		self.obs_data = pd.read_csv(obs_path)
		#read AS test error
		self.dps_Ierr = []
		self.dps_Verr = []
		self.completeness_V = []
		self.completeness_I = []
		for i in range(80):
			dp_Ierr = pd.read_csv("{}/Ierr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=3,names=['Ierr'])
			dp_Verr = pd.read_csv("{}/Verr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=3,names=['Verr'])
			i_err, cdf = self.empirical_cdf(dp_Ierr['Ierr'].values)
			self.dps_Ierr.append(interp1d(cdf, i_err, bounds_error=False, fill_value='extrapolate'))
			v_err, cdf = self.empirical_cdf(dp_Verr['Verr'].values)
			self.dps_Verr.append(interp1d(cdf, v_err, bounds_error=False, fill_value='extrapolate'))
			self.completeness_V.append(pd.read_csv("{}/Verr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=1,names=["#","Npts","Radius","Mag","Completeness"],nrows=1)['Completeness'].values[0])
			self.completeness_I.append(pd.read_csv("{}/Ierr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=1,names=["#","Npts","Radius","Mag","Completeness"],nrows=1)['Completeness'].values[0])
		df_V_bound = pd.read_csv("{}/VerrBoundary.dat".format(phot_path),sep='\s+',skiprows=3, names=['Rad', 'Mag', 'Nstar', 'Completness'])
		self.Rad_bounds = np.sort(np.unique(df_V_bound['Rad'].values))
		self.V_mag_bounds = np.sort(np.unique(df_V_bound['Mag'].values))
		df_I_bound = pd.read_csv("{}/IerrBoundary.dat".format(phot_path),sep='\s+',skiprows=3, names=['Rad', 'Mag', 'Nstar', 'Completness'])
		self.Rad_bounds = np.sort(np.unique(df_I_bound['Rad'].values))
		self.I_mag_bounds = np.sort(np.unique(df_I_bound['Mag'].values))
		#read fiducial lines
		fiducial = pd.read_csv(fiducial_path,names=['vi','v','c5','c95','perp'])
		self.V_fid = fiducial['v'].values
		I_fid = fiducial['v'].values - fiducial['vi'].values
		self.ff = interp1d(self.V_fid, I_fid, bounds_error=False, fill_value='extrapolate')

	#maginalized the completeness from the result of AS test. We assume the radius and magnitude are indepedent
	def marginalize(self):
		self.completeness_V_marginalzied = []
		self.completeness_I_marginalzied = []
		self.completeness_R_marginalized = []
		x = self.obs_data['x'] - self.OBS_center[0]
		y = self.obs_data['y'] - self.OBS_center[1]
		r = np.sqrt(x**2 + y**2)
		r, r_cdf = self.empirical_cdf(r)

		#find r bounds
		i = 0
		k = 0
		idxs_r = []
		while k < (len(r) - 1):
			if (r[k+1] > self.Rad_bounds[i]) & (r[k] <= self.Rad_bounds[i]):
				idxs_r.append(k+1)
				i += 1
				if i == (len(self.Rad_bounds) - 1):
					break
				k += 1
			else:
				k += 1
		#in case it runs outside of the bound
		if len(idxs_r) < len(self.Rad_bounds):
			idxs_r.append(-1)

		#do Vmag
		for j in range(len(self.V_mag_bounds)):
			maginalzied_prob = 0
			for i in range(len(self.Rad_bounds)):
				if i==0:
					maginalzied_prob += self.completeness_V[len(self.Rad_bounds)*j + i]*r_cdf[idxs_r[i]]
				else:
					maginalzied_prob += self.completeness_V[len(self.Rad_bounds)*j + i]*(r_cdf[idxs_r[i]] - r_cdf[idxs_r[i-1]])
			self.completeness_V_marginalzied.append(maginalzied_prob)

		#do Imag
		for j in range(len(self.I_mag_bounds)):
			maginalzied_prob = 0
			for i in range(len(self.Rad_bounds)):
				if i==0:
					maginalzied_prob += self.completeness_I[len(self.Rad_bounds)*j + i]*r_cdf[idxs_r[i]]
				else:
					maginalzied_prob += self.completeness_I[len(self.Rad_bounds)*j + i]*(r_cdf[idxs_r[i]] - r_cdf[idxs_r[i-1]])
			self.completeness_I_marginalzied.append(maginalzied_prob)
		return r, r_cdf

	#find the completeness with respect to radius from the AS test
	def find_r_completness(self,as_test_path,pix_miss=0.5,mag_miss=0.75):
		AS_test = pd.read_csv(as_test_path,sep='\s+',names=['x_in','y_in','v_in','i_in','x_out','y_out','v_out','i_out'],skiprows=1)
		r_in = np.sqrt((AS_test['x_in'] - self.AS_center[0])**2 + (AS_test['y_in'] - self.AS_center[1])**2)
		r_out = np.sqrt((AS_test['x_out'] - self.AS_center[0])**2 + (AS_test['y_out'] - self.AS_center[1])**2)
		r_diff = np.abs(r_in - r_out)
		v_diff = np.abs(AS_test['v_in'] - AS_test['v_out'])
		i_diff = np.abs(AS_test['i_in'] - AS_test['i_out'])
		input_num = np.zeros(len(self.Rad_bounds))
		good_num = np.zeros(len(self.Rad_bounds))
		for i in range(len(AS_test)):
			j = 0
			while j < (len(self.Rad_bounds)):
				if r_in[i] < self.Rad_bounds[j]:
					input_num[j] += 1
					if (r_diff[i] < pix_miss) & (v_diff[i] < mag_miss) & (i_diff[i] < mag_miss):
						good_num[j] += 1
					j = np.inf
				else:
					j += 1
		return good_num/input_num
	 
	#inferencing the real cdf from the observational data with completeness test
	def cdf_inference_v(self,v):
		v, cdf = self.empirical_cdf(v)
		len_v = len(v)
		i = self.ff(v)
		cdf_new = np.zeros(len_v)
		for j in range(len_v):
			_, idx = self.find_fit_val(self.V_mag_bounds,v[j])
			v_comp = self.completeness_V_marginalzied[idx]
			_, idx = self.find_fit_val(self.I_mag_bounds,i[j])
			i_comp = self.completeness_I_marginalzied[idx]
			if j == 0:
				cdf_new[0] = cdf[0]/(v_comp*i_comp)
			else:
				cdf_new[j] = cdf_new[j-1] + (cdf[j] - cdf[j-1])/(v_comp*i_comp)
		cdf_new /= cdf_new[-1]
		return v, cdf_new
	
	def get_cdf(self,cdf_path,as_test_path):
		r, r_cdf = self.marginalize()
		v, v_cdf_new = self.cdf_inference_v(self.obs_data['v'].values)
		r_completeness = self.find_r_completness(as_test_path)
		len_r = len(r)
		r_cdf_new = np.zeros(len_r)
		for i in range(len_r):
			_, idx = self.find_fit_val(self.Rad_bounds,r[i])
			r_comp = r_completeness[idx]
			if i == 0:
				r_cdf_new[0] = r_cdf[0]/r_comp
			else:
				r_cdf_new[i] = r_cdf_new[i-1] + (r_cdf[i] - r_cdf[i-1])/r_comp
		r_cdf_new /= r_cdf_new[-1]
		#write files
		d_r = {'cdf': r_cdf_new, 'r': r}
		d_v = {'cdf': v_cdf_new, 'v': v}
		df = pd.DataFrame(data = d_r)
		df.to_csv("{}/marginal_r.csv".format(cdf_path),index=False)
		df = pd.DataFrame(data = d_v)
		df.to_csv("{}/marginal_v.csv".format(cdf_path),index=False)

	#read in marginal files
	def read_marginals(self, cdf_path):
		df_v = pd.read_csv("{}/marginal_v.csv".format(cdf_path))
		self.v_cdf = interp1d(df_v['cdf'].values, df_v['v'].values, bounds_error=False, fill_value='extrapolate')
		df_r = pd.read_csv("{}/marginal_r.csv".format(cdf_path))
		self.r_cdf = interp1d(df_r['cdf'].values, df_r['r'].values, bounds_error=False, fill_value='extrapolate')

	def resample(self, Binary_Fraction, sample_pt, obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,k):
		#sample v magitude and calculate i from the fiducial isochrone
		v_iso = self.v_cdf(np.random.rand(sample_pt))
		i_iso = self.ff(v_iso)
		#convert the magnitude to AS_test magnitude
		v_iso -= self.V_diff
		i_iso -= self.I_diff
		v_min = obs_v_min - self.V_diff
		v_max = obs_v_max - self.V_diff
		vi_min = obs_vi_min - self.V_diff + self.I_diff
		vi_max = obs_vi_max - self.V_diff + self.I_diff
		#sample radius
		r_iso = self.r_cdf(np.random.rand(sample_pt))
		#make binaries (still working out how to make a binary)
		Binary_num = round(sample_pt*Binary_Fraction)

		#prepare the random number we need for completeness test and error
		completeness_test_v = np.random.rand(sample_pt)
		completeness_test_i = np.random.rand(sample_pt)
		err_cdf_v = np.random.rand(sample_pt)
		err_cdf_i = np.random.rand(sample_pt)
		#define good_v and good_i
		good_v = []
		good_vi = []
		#shift sample points based on their corresponding completeness or error
		num_v_bin = len(self.V_mag_bounds)
		num_r_bin = len(self.Rad_bounds)
		num_i_bin = num_v_bin
		for i in range(sample_pt):
			_, v_bin = self.find_fit_val(self.V_mag_bounds,v_iso[i])
			_, r_bin = self.find_fit_val(self.Rad_bounds,r_iso[i])
			if completeness_test_v[i] < self.completeness_V[v_bin*num_r_bin + r_bin]:
				_, i_bin =  self.find_fit_val(self.I_mag_bounds,i_iso[i])
				#completeness test
				if completeness_test_i[i] < self.completeness_I[i_bin*num_r_bin + r_bin]:
					verr = self.dps_Verr[v_bin*num_r_bin + r_bin](err_cdf_v[i])
					ierr = self.dps_Ierr[i_bin*num_r_bin + r_bin](err_cdf_i[i])
					temp_v = v_iso[i] + verr
					temp_i = i_iso[i] + ierr
					temp_vi = temp_v - temp_i
					#data cut test
					if (np.abs(verr) < 0.08) & (temp_v > v_min) & (temp_v < v_max) & (temp_vi > vi_min) & (temp_vi < vi_max):
						good_v.append(temp_v)
						good_vi.append(temp_vi)
		#convert back to observational magnitude
		good_v = np.array(good_v) + self.V_diff
		good_vi = np.array(good_vi) + self.V_diff - self.I_diff
		df = pd.DataFrame(data={'v':good_v, 'vi':good_vi})
		print("Done generating sCMD for {}, get {} of points".format(k,len(df)))
		return df




	def __init__(self, GC_name, start, end, Binary_Fraction = 0.02, write_vorbin=False, Tb_size=30,write_cmd=False, sample_pt=4000000):
		#define boundaris
		if GC_name == 'M55':
			obs_vi_max = 0.792
			obs_vi_min = 0.462
			obs_v_max = 19.28
			obs_v_min = 15.296
		width_coeff = (obs_v_max - obs_v_min)/(obs_vi_max - obs_vi_min)
		#correct the difference between obs data and as test
		if GC_name == 'M55':
			self.AS_center = [3005.49976,2997.75391]
			self.OBS_center = [2995.02393,3020.84351]
			self.V_diff=30.98
			self.I_diff = self.V_diff-0.737
		#define all the path for read and write
		if GC_name == 'M55':
			resample_path = '/media/sf_share/M55_data/resample'
			data_path = '/home/mying/Desktop/GC_Ages/Age-of-GCs/M55_data'
			photometry_folder = '/home/mying/Desktop/GC_Ages/Age-of-GCs/Photometry'
		photometry_path = "{}/M55_inputfiles".format(photometry_folder)
		vorbin_path = "{}/vorbin".format(resample_path )
		chi2_path = "{}/outchi2".format(resample_path)
		cmd_path = "{}/cmd".format(resample_path )
		obs_path = "{}/M55_fitstars.dat".format(data_path)
		as_test_path = "{}/M55artstars.dat".format(data_path)
		fiducial_path = "{}/fiducial_lines.csv".format(data_path)
		cdf_path = "{}/cdf".format(resample_path)
		#check those directories exist
		self.check_file(photometry_path)
		self.check_file(fiducial_path)
		self.check_file(as_test_path)
		self.check_directories(vorbin_path)
		self.check_directories(chi2_path)
		self.check_directories(cmd_path)
		self.check_directories(cdf_path)
		#assign other global variables
		self.Tb_size = Tb_size
		self.sample_pt = sample_pt
		#find both cuts from observational data
		if GC_name == 'M55':
			self.MSTO_cut = 17.0
			self.GB_cut = 18.0
		#read obs data and photometry data
		self.read_input(obs_path,photometry_path,fiducial_path)
		#check if cdf files exits
		if os.path.exists("{}/marginal_r.csv".format(cdf_path)) == False:
			self.get_cdf(cdf_path,as_test_path)
			print('Done writing cdf')
		#read cdf files
		self.read_marginals(cdf_path)
		#run resample
		self.chi2 = []
		for k in range(start, end):
			print("Starting {}th resample".format(k))
			df = self.resample(Binary_Fraction, sample_pt, obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,k)
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
		self.writeout_resample(chi2_path,start,end)
