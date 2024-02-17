import vorbin_chi2
import subprocess
import os
import time
from global_var import sCMD_vars
from path_config import resample_path

class resample_iso:
	def CMD_gen(self,resample_num):
		subprocess.run(['./TestCMDPAR.sh', str(self.mc_num), str(self.mc_num), str(self.pdmf), str(self.binary), '4000000', str(self.age), str(self.GC_name),str(self.feh), '2.0'])
		source = self.outcmd_path + "/mc{}.a{}".format(self.mc_num,self.age)
		dest = self.outcmd_path + "/mc{}.a{}_{}".format(self.mc_num,self.age, str(resample_num))
		os.rename(source,dest)

	def CMD_check(self,resample_num):
		check_time = 0
		exist = False
		while check_time < 10 and exist == False:
			check_time += 1
			if os.path.exists("{}/mc{}.a{}_{}".format(self.outcmd_path,self.mc_num, self.age,resample_num)):
				exist = True
				break
			else:
				#sleep for 10 second
				time.sleep(10)
		if exist == False:
			self.CMD_gen(resample_num)
	
	def run_vorbin(self, obs_i, resample_i):
		vorbin_chi2.chi2_iso(self.GC_name, self.mc_num, self.age, obs_i, resample_i, UniSN=True, write_vorbin=True)

	def rmCMD(self):
		#remove searched sCMD file
		file = "mc{}.a{}".format(self.mc_num, self.age)
		path = os.path.join(self.outcmd_path, file)
		os.remove(path)

	def __init__(self, GC_name, mc_num, age, resample_num=100):
		self.GC_name = str(GC_name)
		#define global variables
		self.outcmd_path ="{}{}/outcmd".format(resample_path, self.GC_name)
		#define globar variables used to generate sCMD
		self.feh, self.binary, self.pdmf = sCMD_vars(GC_name)
		self.mc_num = str(mc_num)
		if float(age) < 10000:
			self.age = '0'+str(age)
		else:
			self.age = str(age)
		#time it
		start_time = time.time()
		#generate CMD
		for i in range(resample_num):
			self.CMD_gen(i)
			self.CMD_check(i)
		#run analysis
		for i in range(resample_num):
			for j in range(resample_num):
				if i != j:
					self.run_vorbin(i,j)
		# self.rmCMD()
		print("Done MC{} Age{}".format(self.mc_num, self.age))
		print("--- %s seconds ---" % (time.time() - start_time))