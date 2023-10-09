#This code requires three inputs: GC_name, mc_num and age
#This code outputs a list of chi2 for dm and reddening combination
import numpy as np
import pandas as pd
import os
from vorbin.voronoi_2d_binning import voronoi_2d_binning
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.optimize import differential_evolution as DE
from skopt import gp_minimize
from scipy.interpolate import LinearNDInterpolator
#from bayes_opt import BayesianOptimization
#from skopt import gp_minimize
#from skopt.space.space import Real

#force the code to run on 1core as serial jobs on Stampede2
#os.environ['OPENBLAS_NUM_THREADS'] = '1'
#os.environ['NUMEXPR_NUM_THREADS'] = '1'
#os.environ['MKL_NUM_THREADS'] = '1'

class utiles:
	#Divide-and-Conquer method to find the 2d ecdf
	def _rank2(self, points, mask=None):
		N = points.shape[0]
		N2 = N//2
		if N == 1:
			return 0
		else:
			idx = np.argpartition(points[:,0], N2)
			idxA_ = idx[:N2]
			idxA = np.zeros(N, dtype=bool)
			idxA[idxA_] = True
			if mask is not None:
				NAm = np.sum(idxA & mask)
				points_reduced = np.vstack((points[idxA & mask], points[~idxA & ~mask]))
			else:
				NAm = np.sum(idxA)
				points_reduced = np.vstack((points[idxA], points[~idxA]))
			count_points = np.zeros(points_reduced.shape[0], dtype=bool)
			count_points[:NAm] = True
			idxY = np.argsort(points_reduced[:,1])
			idxYr = np.zeros_like(idxY)
			idxYr[idxY] = np.arange(idxY.shape[0]) # inverse of idxY
			count_points = count_points[idxY]
			numA = np.cumsum(count_points)[idxYr]
			rank = np.zeros(N, dtype=int)
			if mask is not None:
				rank[idxA] = self._rank2(points[idxA], mask[idxA])
				rank[~idxA] = self._rank2(points[~idxA], mask[~idxA])
				rank[~idxA & ~mask] += numA[NAm:]
			else:
				rank[idxA] = self._rank2(points[idxA])
				rank[~idxA] = self._rank2(points[~idxA])
				rank[~idxA] += numA[NAm:]
			return rank

	def rankn(self, points, mask=None):
		N = points.shape[0]
		N2 = N//2
		if mask is None:
			mask = np.ones(N, dtype=bool)
			first_call = True
		else:
			first_call = False
		if N == 1:
			return 0
		if points.shape[1] == 2:
			if first_call:
				return self._rank2(points)
			else:
				return self._rank2(points, mask)
		idx = np.argpartition(points[:,0], N2)
		idxA_ = idx[:N2]
		idxA = np.zeros(N, dtype=bool)
		idxA[idxA_] = True
		rank = np.zeros(N, dtype=int)
		rank[idxA] = self.rankn(points[idxA], mask[idxA])
		rank[~idxA] = self.rankn(points[~idxA], mask[~idxA]) + self.rankn(points[:,1:], idxA*mask)[~idxA]
		return rank
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
		#search through observed data
		#Utilize Bayesian optimization method to find the best fit DM and red value
		# Bounded region of parameter space
		# pbounds = {'dm': (dm_min, dm_max), 'red': (red_min, red_max)}

		# optimizer = BayesianOptimization(
		# 	f=self.dm_red_search,
		# 	pbounds=pbounds,
		# 	random_state=1,
		# 	verbose=0,
		# )

		# if write_chi2_log == True:
		# 	from bayes_opt.logger import JSONLogger
		# 	from bayes_opt.event import Events
		# 	logger = JSONLogger(path="{}/mc{}_age{}_logs.log".format(chi2_path,self.mc_num,self.iso_age))
		# 	optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
		
		# optimizer.maximize(
		# 	init_points=100,
		# 	n_iter=100,
		# )
		# chi2 = np.array([self.iso_age, optimizer.max['params']['dm'], optimizer.max['params']['red'] ,-optimizer.max['target'], len(self.obs_cut(optimizer.max['params']['dm'], optimizer.max['params']['red']))])
		
		#Try to use Bayesian optimization from skopt
		# dim1 = Real(name='dm', low=dm_min, high=dm_max)
		# dim2 = Real(name='red', low=red_min, high=red_max)
		# bounds = [dim1, dim2]
		# res = gp_minimize(self.dm_red_search,                  # the function to minimize
        #           bounds,      # the bounds on each dimension of x
        #           acq_func="EI",
        #           n_calls=100,         # the number of evaluations of f
        #           n_initial_points=10,
        #           noise=0,
        #           acq_optimizer="lbfgs",
        #           n_restarts_optimizer=5,
        #           #n_random_starts=100,  # the number of random initialization points
        #          )
		# dm_fit = res['x'][0]
		# red_fit = res['x'][1]
		# chi2_fit = res['fun']
		# retval = np.array([self.iso_age, dm_fit, red_fit, chi2_fit, len(self.obs_cut(dm_fit, red_fit))])


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
				chi2 = self.dm_red_search([dm,red])
				if chi2 < chi2_fit:
					chi2_fit = chi2
					dm_fit = dm
					red_fit = red

		#write out the result
		retval = np.array([self.iso_age, dm_fit, red_fit, chi2_fit, len(self.obs_cut(dm_fit, red_fit))])
		#print(chi2)
		pd.DataFrame(retval).to_csv("{}/chi2_a{}_mc{}".format(chi2_path,self.iso_age,self.mc_num),header=None, index=None)
		# for dm in dms:
		# 	for red in reds:
		# 		obs_cut = self.obs_data[(self.obs_data['vi'] - red < (obs_vi_max - reds[-1])) & (self.obs_data['vi'] - red > (obs_vi_min - reds[0]))& (self.obs_data['v'] - dm < (obs_v_max - dms[-1])) & (self.obs_data['v'] - dm > (obs_v_min - dms[0]))]
		# 		obs_size = len(obs_cut)
		# 		vi_32 = np.float32((obs_cut['vi'].values - red)*width_coeff)
		# 		v_32 = np.float32(obs_cut['v'].values - dm)
		# 		bin_count = self.search_vorbin(XBar, YBar*width_coeff, obs_size, vi_32, v_32)
		# 		#calculate chi2
		# 		chi2.append([age, dm, red, np.inner(np.divide(bin_count,bin_count_std/(total_pt/obs_size)) - 1, bin_count - bin_count_std/(total_pt/obs_size))])
		# self.chi2 = chi2

	def __init__(self, GC_name, mc_num, iso_age, UniSN=False, write_vorbin=False, Tb_size=30):
		#define distance modulus and reddening ranges
		if GC_name == 'M55':
			# dm_max = 14.40
			# dm_min = 13.40
			# red_max = 0.20
			# red_min = 0.00
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
		self.main(write_vorbin,cmd_path, dm_max, dm_min, red_max, red_min,iso_path,chi2_path,UniSN=UniSN)
		print("done mc{}".format(self.mc_num))

class KS_2d(utiles):
	def read_input(self,path):
		#read M92 observed data
		obs_data = pd.read_csv(path)
		# names = ['v','v_err','i','i_err','vi','vi_err','x','y']
		# for i in range(len(names)):
		# 	for j in range(len(obs_data.columns)):
		# 		if names[i] == obs_data.columns[j]:
		# 			setattr(self,names[i] + '_idx', j)
		self.obs_data = pd.read_csv(path)[['vi','v']].to_numpy()
		self.obs_data[:,1] = - self.obs_data[:,1]

	def evaluate(self, theta):
		dm,red = theta
		#cut the sCMD data first
		Obs_DM_Red_CR = self.obs_data - [red, -dm]
		sCMD_in_range = self.sCMD[(self.sCMD['v'] < self.obs_v_max - dm) & (self.sCMD['v'] > self.obs_v_min - dm)]
		fit_values = sCMD_in_range.values[:self.num_stars]
		fit_values[:,1] = -fit_values[:,1]
		sCMD_emp_ff = LinearNDInterpolator(fit_values, self.rankn(fit_values)/len(fit_values), fill_value=0)
		sCMD_Pred_cdf = sCMD_emp_ff(Obs_DM_Red_CR)
		mask = (self.Obs_cdf != 0.0) & (sCMD_Pred_cdf != 0.0)
		return np.max(np.abs(sCMD_Pred_cdf[mask] - self.Obs_cdf[mask]))

	def main(self, cmd_path, dm_max, dm_min, red_max, red_min,chi2_path):
		age = self.iso_age
		#read cmd files
		self.sCMD = pd.read_csv("{}/mc{}.a{}".format(cmd_path,self.mc_num,age),sep='\s+',names=['vi','v'],skiprows=3)
		#go through the search process
		self.obs_v_max = max(-self.obs_data[:,1])
		self.obs_v_min = min(-self.obs_data[:,1])
		self.Obs_cdf = self.rankn(self.obs_data)/len(self.obs_data)
		res = gp_minimize(self.evaluate,                  # the function to minimize
                  [(dm_min, dm_max), (red_min, red_max)],      # the bounds on each dimension of x
                  acq_func="EI",      # the acquisition function
                  n_calls=50,         # the number of evaluations of f
                  n_random_starts=10,
                  verbose=False)
		retval = np.array([self.iso_age, res.x[0], res.x[1], res.fun])
		#print(chi2)
		pd.DataFrame(retval).to_csv("{}/chi2_a{}_mc{}".format(chi2_path,self.iso_age,self.mc_num),header=None, index=None)

	#use 2d KS test to calculate the Metric for the input isochrones
	def __init__(self, GC_name, mc_num, iso_age,num_stars=60000):
		#define distance modulus and reddening ranges
		if GC_name == 'M55':
			dm_max = 14.1
			dm_min = 13.8
			red_max = 0.15
			red_min = 0.08
		#define other global variables
		self.num_stars=num_stars
		self.mc_num = str(mc_num)
		self.iso_age = str(iso_age)
		if GC_name == 'M92':
			self.feh = 230
		elif GC_name == 'M55':
			self.feh = 190
		#define all the path for read and write
		obs_data_path = "/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/{}/simulateCMD/{}_fitstars_ZPCR.dat".format(GC_name,GC_name)
		vorbin_path = "/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/{}/vorbin".format(GC_name)
		self.vorbin_path = vorbin_path
		chi2_path = "/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/{}/outchi2".format(GC_name)
		cmd_path = "/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/{}/simulateCMD/outcmd".format(GC_name)
		iso_path = "/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/{}/outiso".format(GC_name)
		#check those directories exist
		self.check_file(obs_data_path)
		self.check_directories(vorbin_path)
		self.check_directories(chi2_path)
		self.check_directories(cmd_path)
		self.check_directories(iso_path)
		#run code
		self.read_input(obs_data_path)
		self.main(cmd_path, dm_max, dm_min, red_max, red_min,chi2_path)
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
			self.dps_Ierr.append(pd.read_csv("{}/Ierr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=3,names=['Ierr']))
			self.dps_Verr.append(pd.read_csv("{}/Verr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=3,names=['Verr']))
			self.completeness_V.append(pd.read_csv("{}/Verr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=1,names=["#","Npts","Radius","Mag","Completeness"],nrows=1)['Completeness'].values[0])
			self.completeness_I.append(pd.read_csv("{}/Ierr{:02d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=1,names=["#","Npts","Radius","Mag","Completeness"],nrows=1)['Completeness'].values[0])

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
			path = "{}/resample_{}".format(path,k)
			df_resample.to_csv(path,index=False)
		return df_resample


	def __init__(self, GC_name, start, end, MSTO_cut, GB_cut, write_vorbin=False, Tb_size=30,write_cmd=False, sample_pt=4000000):
		#define boundaris
		if GC_name == 'M55':
			obs_vi_max = 0.792
			obs_vi_min = 0.462
			obs_v_max = 19.28
			obs_v_min = 15.296
		width_coeff = (obs_v_max - obs_v_min)/(obs_vi_max - obs_vi_min)
		#define all the path for read and write
		if GC_name == 'M55':
			# repo_path = '/home/mying/Desktop/GC_Ages/Age-of-GCs'
			# resample_path = '/media/sf_share/{}_data/resample'.format(GC_name)
			repo_path = '/home/mying/Desktop/Repositories/Age-of-GCs'
			resample_path = "/home/mying/Desktop/{}_data/resample".format(GC_name)
		data_path = "{}/{}_data".format(repo_path, GC_name)
		photometry_folder = "{}/Photometry".format(repo_path)
		photometry_path = "{}/{}_inputfiles".format(photometry_folder,GC_name)
		vorbin_path = "{}/vorbin".format(resample_path )
		chi2_path = "{}/outchi2".format(resample_path)
		cmd_path = "{}/cmd".format(resample_path )
		obs_data_path = "{}/{}_fitstars_with_bins.dat".format(data_path,GC_name)
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
		chi2 = []
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
			chi2.append([k, np.inner(np.divide(bin_count,bin_count_std/(total_pt/obs_size)) - 1, bin_count - bin_count_std/(total_pt/obs_size))])  
		self.writeout_resample(chi2_path,start,end,chi2)

#this is similar to resample class but utilizing the fiducial isochrone generated from fidanka
class resample_fidanka(utiles):
	def read_input(self,obs_path,phot_path,fiducial_path,binary_path,obs_vi_max,obs_vi_min,obs_v_max,obs_v_min):
		#read observed data
		self.obs_data = pd.read_csv(obs_path)
		self.obs_data = self.obs_data[(self.obs_data['v'] <= obs_v_max) & (self.obs_data['v'] >= obs_v_min) & (self.obs_data['vi'] <= obs_vi_max) & (self.obs_data['vi'] >= obs_vi_min)]
		#read AS test error
		self.dps_Ierr = []
		self.dps_Verr = []
		self.completeness_V = []
		self.completeness_I = []
		for i in range(140):
			dp_Ierr = pd.read_csv("{}/Ierr{:03d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=3,names=['Ierr'])
			dp_Verr = pd.read_csv("{}/Verr{:03d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=3,names=['Verr'])
			i_err, cdf = self.empirical_cdf(dp_Ierr['Ierr'].values)
			self.dps_Ierr.append(interp1d(cdf, i_err, bounds_error=False, fill_value='extrapolate'))
			v_err, cdf = self.empirical_cdf(dp_Verr['Verr'].values)
			self.dps_Verr.append(interp1d(cdf, v_err, bounds_error=False, fill_value='extrapolate'))
			self.completeness_V.append(pd.read_csv("{}/Verr{:03d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=1,names=["#","Npts","Radius","Mag","Completeness"],nrows=1)['Completeness'].values[0])
			self.completeness_I.append(pd.read_csv("{}/Ierr{:03d}s.dat".format(phot_path,i + 1),sep='\s+',skiprows=1,names=["#","Npts","Radius","Mag","Completeness"],nrows=1)['Completeness'].values[0])
		self.completeness_V = np.array(self.completeness_V)
		self.completeness_I = np.array(self.completeness_I)
		self.dps_Verr = np.array(self.dps_Verr)
		self.dps_Ierr = np.array(self.dps_Ierr)
		df_V_bound = pd.read_csv("{}/VerrBoundary.dat".format(phot_path),sep='\s+',skiprows=3, names=['Rad', 'Mag', 'Nstar', 'Completness'])
		self.Rad_bounds = np.sort(np.unique(df_V_bound['Rad'].values))
		self.V_mag_bounds = np.sort(np.unique(df_V_bound['Mag'].values))
		df_I_bound = pd.read_csv("{}/IerrBoundary.dat".format(phot_path),sep='\s+',skiprows=3, names=['Rad', 'Mag', 'Nstar', 'Completness'])
		self.Rad_bounds = np.sort(np.unique(df_I_bound['Rad'].values))
		self.I_mag_bounds = np.sort(np.unique(df_I_bound['Mag'].values))
		#read fiducial lines
		fiducial = pd.read_csv(fiducial_path,dtype=float)
		self.V_fid = fiducial['v'].values
		I_fid = fiducial['v'].values - fiducial['vi'].values
		self.ff = interp1d(self.V_fid, I_fid, bounds_error=False, fill_value='extrapolate')
		#read binary charts
		binary_data = pd.read_csv(binary_path)
		v_len = 600
		q_len = 100
		v = np.linspace(0,6,v_len)
		q = np.linspace(0.5,1,q_len)
		v_diff = binary_data['v_diff'].values.reshape(v_len,q_len)
		vi_diff = binary_data['vi_diff'].values.reshape(v_len,q_len)
		self.f_vdiff = RectBivariateSpline(v, q, v_diff)
		self.f_vidiff = RectBivariateSpline(v, q, vi_diff)


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

	def resample(self, Binary_Fraction, sample_pt, obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,k,cmd_path,write_cmd):
		#sample v magitude and calculate i from the fiducial isochrone
		v_iso = self.v_cdf(np.random.rand(sample_pt))
		i_iso = self.ff(v_iso)
		vi_iso = v_iso - i_iso
		#sample radius
		r_iso = self.r_cdf(np.random.rand(sample_pt))
		#test generating binary from random sampling the fiducial isochrone
		#make binaries (still working out how to make a binary)
		Binary_num = round(sample_pt*Binary_Fraction)
		q_iso = 0.5 + np.random.rand(Binary_num)/2
		#find the change 
		v_diff = np.zeros(sample_pt)
		i_diff = np.zeros(sample_pt)
		vi_diff = np.zeros(sample_pt)
		v_diff[:Binary_num] = self.f_vdiff(v_iso[:Binary_num] - (self.V_SGB - self.mag_cut), q_iso,grid=False)
		vi_diff[:Binary_num] = self.f_vidiff(v_iso[:Binary_num] - (self.V_SGB - self.mag_cut), q_iso,grid=False)
		i_diff =  v_diff - vi_diff
		v_iso[:Binary_num] = v_iso[:Binary_num] + v_diff[:Binary_num]
		vi_iso[:Binary_num] = vi_iso[:Binary_num] + vi_diff[:Binary_num]
		i_iso[:Binary_num] = v_iso[:Binary_num] - vi_iso[:Binary_num]
		# v_binary = []
		# i_binary = []
		# for i in range(Binary_num):
		# 	v_sec = self.v_cdf(np.random.rand())
		# 	while v_sec < v_iso[i]:
		# 		v_sec = self.v_cdf(np.random.rand())
		# 	v_binary.append(v_sec)
		# 	i_binary.append(self.ff(v_sec))
		# v_binary = np.array(v_binary)
		# i_binary = np.array(i_binary)
		# # v_binary = self.v_cdf(np.random.rand(Binary_num))
		# # i_binary = self.ff(v_binary)
		# # v_binary -= self.V_diff
		# # i_binary -= self.I_diff
		# v_iso[:Binary_num] = -2.5*np.log10(np.power(10,-0.4*v_iso[:Binary_num]) + np.power(10,-0.4*v_binary))
		# i_iso[:Binary_num] = -2.5*np.log10(np.power(10,-0.4*i_iso[:Binary_num]) + np.power(10,-0.4*i_binary))
		#convert the magnitude to AS_test magnitude
		v_iso -= self.V_diff
		i_iso -= self.I_diff
		v_min = obs_v_min - self.V_diff
		v_max = obs_v_max - self.V_diff
		vi_min = obs_vi_min - self.V_diff + self.I_diff
		vi_max = obs_vi_max - self.V_diff + self.I_diff
		#prepare the random number we need for completeness test and error
		completeness_test_v = np.random.rand(sample_pt)
		completeness_test_i = np.random.rand(sample_pt)
		err_cdf_v = np.random.rand(sample_pt)
		err_cdf_i = np.random.rand(sample_pt)
		#shift sample points based on their corresponding completeness or error
		num_v_bin = len(self.V_mag_bounds)
		num_r_bin = len(self.Rad_bounds)
		v_bin = np.clip(np.searchsorted(self.V_mag_bounds, v_iso, side='left'),0,len(self.V_mag_bounds) - 1)
		r_bin = np.clip(np.searchsorted(self.Rad_bounds, r_iso, side='left'),0,len(self.Rad_bounds) - 1)
		i_bin = np.clip(np.searchsorted(self.I_mag_bounds, i_iso, side='left'),0,len(self.I_mag_bounds) - 1)
		#completenes test results
		v_bin_full = v_bin*num_r_bin + r_bin
		i_bin_full = i_bin*num_r_bin + r_bin
		v_result = completeness_test_v < self.completeness_V[v_bin_full]
		i_result = completeness_test_i < self.completeness_I[i_bin_full]
		#shift the point by error from AS test
		verr = np.zeros(sample_pt)
		ierr = np.zeros(sample_pt)
		for i in range(len(self.dps_Verr)):
			mask = v_bin_full == i
			verr[mask] = self.dps_Verr[i](err_cdf_v[mask])
			mask = i_bin_full == i
			ierr[mask] = self.dps_Ierr[i](err_cdf_i[mask])
		temp_v = v_iso + verr
		temp_i = i_iso + ierr
		temp_vi = temp_v - temp_i
		#combine all conditions
		All_tests = (np.abs(ierr+i_diff) < 0.08) & (np.abs(verr+v_diff) < 0.08) & (temp_v > v_min) & (temp_v < v_max) & (temp_vi > vi_min) & (temp_vi < vi_max) & v_result & i_result
		#select points satisfied all conditions
		good_v = temp_v[All_tests]
		good_vi = temp_vi[All_tests]
		# for i in range(sample_pt):
		# 	if v_result[i] & i_result[i]:
		# 		verr = self.dps_Verr[v_bin*num_r_bin + r_bin](err_cdf_v[i])
		# 		ierr = self.dps_Ierr[i_bin*num_r_bin + r_bin](err_cdf_i[i])
		# 		temp_v = v_iso[i] + verr
		# 		temp_i = i_iso[i] + ierr
		# 		temp_vi = temp_v - temp_i
		# 		#data cut test
		# 		if (np.abs(verr) < 0.08) & (temp_v > v_min) & (temp_v < v_max) & (temp_vi > vi_min) & (temp_vi < vi_max):
		# 			good_v.append(temp_v)
		# 			good_vi.append(temp_vi)
		#convert back to observational magnitude
		good_v = np.array(good_v) + self.V_diff
		good_vi = np.array(good_vi) + self.V_diff - self.I_diff
		df = pd.DataFrame(data={'v':good_v, 'vi':good_vi})
		print("Done generating sCMD for {}, get {} of points".format(k,len(df)))
		if write_cmd == True:
			path = "{}/resample_{}".format(cmd_path,k)
			df.to_csv(path,index=False)
		return df

	def calculate_chi2(self,Binary_Fraction, sample_pt, obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,vorbin_path,chi2_path,start, end, write_vorbin,cmd_path,write_cmd):
		chi2 = []
		for k in range(start, end):
			print("Starting {}th resample".format(k))
			df = self.resample(Binary_Fraction, sample_pt, obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,k,cmd_path,write_cmd)
			total_pt = len(df)
			V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB = self.divide_data(df)
			x_gen, y_gen = self.generate_vorbin([V_MS, VI_MS, V_MSTO, VI_MSTO, V_GB, VI_GB])
			#reduce memory usage for matrix operations
			XBar = np.float32(x_gen)
			YBar = np.float32(y_gen)
			v_32 = np.float32(df['v'].values)
			vi_32 = np.float32(df['vi'].values*self.width_coeff)
			#find standard bin count by search through all the theoretical data points
			bin_count_std = self.search_vorbin(XBar, YBar, total_pt, vi_32, v_32)
			if write_vorbin == True:
				#no age for resample
				age = 0
				self.writevorbin(x_gen, y_gen, bin_count_std, k, age, vorbin_path)
			#fit observation data
			obs_size = len(self.obs_data)
			vi_32 = np.float32(self.obs_data['vi'].values*self.width_coeff)
			v_32 = np.float32(self.obs_data['v'].values)
			bin_count = self.search_vorbin(XBar, YBar, obs_size, vi_32, v_32)
			#calculate chi2
			chi2.append([k, np.inner(np.divide(bin_count,bin_count_std/(total_pt/obs_size)) - 1, bin_count - bin_count_std/(total_pt/obs_size))/len(self.obs_data)])  
		self.writeout_resample(chi2_path,start,end,chi2)



	def __init__(self, GC_name, start, end, write_vorbin=False, Tb_size=30,write_cmd=False, sample_pt=4000000, pool= False):
		#define boundaris
		if GC_name == 'M55':
			self.V_SGB = 17.28
			self.mag_cut = 3
			Binary_Fraction = 0.04
			obs_vi_max = 0.792
			obs_vi_min = 0.462
			obs_v_max = 19.28
			obs_v_min = 15.296
		self.width_coeff = (obs_v_max - obs_v_min)/(obs_vi_max - obs_vi_min)
		#correct the difference between obs data and as test
		if GC_name == 'M55':
			self.AS_center = [3005.49976,2997.75391]
			self.OBS_center = [2995.02393,3020.84351]
			self.V_diff=30.98
			self.I_diff = self.V_diff-0.737
		#define all the path for read and write
		if GC_name == 'M55':
			# repo_path = '/home/mying/Desktop/GC_Ages/Age-of-GCs'
			# resample_path = '/media/sf_share/{}_data/resample'.format(GC_name)
			repo_path = '/home/mying/Desktop/GC_Ages/Age-of-GCs'
			resample_path = "/home/mying/Desktop/ipynb/{}_data/resample".format(GC_name)
		data_path = "{}/{}_data".format(repo_path, GC_name)
		binary_path = "{}/{}_binary_chart.csv".format(data_path, GC_name)
		photometry_folder = "{}/Photometry".format(repo_path)
		photometry_path = "{}/{}_inputfiles".format(photometry_folder,GC_name)
		vorbin_path = "{}/vorbin".format(resample_path )
		chi2_path = "{}/outchi2".format(resample_path)
		cmd_path = "{}/cmd".format(resample_path )
		obs_path = "{}/{}_fitstars.dat".format(data_path,GC_name)
		as_test_path = "{}/{}artstars.dat".format(data_path,GC_name)
		fiducial_path = "{}/fiducial_lines.csv".format(data_path)
		cdf_path = "{}/cdf".format(resample_path)
		#check those directories exist
		self.check_file(binary_path)
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
		self.read_input(obs_path,photometry_path,fiducial_path,binary_path,obs_vi_max,obs_vi_min,obs_v_max,obs_v_min)
		#check if cdf files exits
		if os.path.exists("{}/marginal_r.csv".format(cdf_path)) == False:
			self.get_cdf(cdf_path,as_test_path)
			print('Done writing cdf')
		#read cdf files
		self.read_marginals(cdf_path)
		#run resample
		if pool == False: 
			self.calculate_chi2(Binary_Fraction, sample_pt, obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,vorbin_path,chi2_path, start, end, write_vorbin,cmd_path,write_cmd)
		else:
			from multiprocessing import Pool
			paramlist = []
			total_resample = end - start
			batch_num = total_resample/10
			i = 0
			while i < batch_num:
				paramlist.append((Binary_Fraction, sample_pt, obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,vorbin_path,chi2_path, start + i*10, start + (i+1)*10, write_vorbin,cmd_path,write_cmd))
				i += 1
			paramlist.append((Binary_Fraction, sample_pt, obs_vi_max,obs_vi_min,obs_v_max,obs_v_min,vorbin_path,chi2_path, start + (i-1)*10, end, write_vorbin,cmd_path,write_cmd))
			with Pool(pool, initializer=np.random.seed) as MP_pool:
				MP_pool.starmap(self.calculate_chi2, paramlist)

class chi2_iso(utiles):

	def read_input(self,path):
		names = ['vi','v']
		obs_data = pd.read_csv(path,skiprows=3, sep='\s+',names=names)
		for i in range(len(names)):
			for j in range(len(obs_data.columns)):
				if names[i] == obs_data.columns[j]:
					setattr(self,names[i] + '_idx', j)
		self.obs_data = pd.read_csv(path,skiprows=3, sep='\s+',names=names).to_numpy()

	def __init__(self, GC_name, mc_num, age, obs_i, resample_i, UniSN=False, write_vorbin=False, Tb_size=30):
		#define distance modulus and reddening ranges
		if GC_name == 'M55':
			# dm_max = 14.40
			# dm_min = 13.40
			# red_max = 0.20
			# red_min = 0.00
			dm_max = 14.1
			dm_min = 13.8
			red_max = 0.15
			red_min = 0.08
		if GC_name == 'M92':
			self.feh = 230
		elif GC_name == 'M55':
			self.feh = 190
		self.Tb_size = Tb_size
		#define all the path for read and write
		obs_data_path = "/home/mying/Desktop/{}_resample/outcmd/mc{}.a{}_{}".format(GC_name, mc_num,age, str(obs_i))
		vorbin_path = "/home/mying/Desktop/{}_resample/vorbin/mc{}.a{}_{}".format(GC_name, mc_num,age, str(obs_i))
		chi2_path = "/home/mying/Desktop/{}_resample/outchi2/mc{}.a{}".format(GC_name, mc_num,age)
		cmd_path = "/home/mying/Desktop/{}_resample/outcmd/mc{}.a{}_{}".format(GC_name, mc_num,age, str(resample_i))
		iso_path = "/home/mying/Desktop/{}_resample/outiso/".format(GC_name)
		#check those directories exist
		self.check_file(obs_data_path)
		#self.check_directories(vorbin_path)
		#self.check_directories(chi2_path)
		#self.check_directories(cmd_path)
		self.check_directories(iso_path)
		self.read_input(obs_data_path)
		#read cmd files
		cmd = pd.read_csv(cmd_path,sep='\s+',names=['vi','v'],skiprows=3)
		#go through the search process
		self.obs_vi_max = max(cmd['vi'].values)
		self.obs_vi_min = min(cmd['vi'].values)
		self.obs_v_max = max(cmd['v'].values)
		self.obs_v_min = min(cmd['v'].values)
		self.width_coeff = (self.obs_v_max - self.obs_v_min)/(self.obs_vi_max - self.obs_vi_min)
		#check whether vorbin already existed or not
		if os.path.isfile(vorbin_path) == True:
			vorbin = pd.read_csv(vorbin_path, skiprows=1,names=['x_gen', 'y_gen','bin_count_std'])
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
			self.v_32 = np.float32(cmd['v'].values)
			self.vi_32 = np.float32(cmd['vi'].values*self.width_coeff)
			#find standard bin count by search through all the theoretical data points
			self.total_pt = len(cmd)
			self.bin_count_std = self.search_vorbin(self.XBar, self.YBar, self.total_pt, self.vi_32, self.v_32)
			#write vorbin infor if desired
			if write_vorbin == True:
				bin_loc = np.vstack((self.XBar,self.YBar,self.bin_count_std)).T
				dp = pd.DataFrame(data=bin_loc, columns = ['xNode', 'yNode','bin_count_std'])
				dp.to_csv(vorbin_path, index=False)
		#calculate chi2
		dm,red = 0.0,0.0
		chi2_fit = self.dm_red_search([dm,red])

		#write out the result
		df_retval = pd.DataFrame({'chi2': [chi2_fit], 'obs_i': [obs_i], 'fit_i':[resample_i]})
		df_retval.to_csv(chi2_path,index=False,mode='a',header=not os.path.exists(chi2_path))

		print("done obs:{}, reasample:{}".format(obs_i, resample_i))