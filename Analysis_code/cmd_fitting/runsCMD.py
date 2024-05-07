import vorbin_chi2
import subprocess
import os
import time
from global_var import sCMD_vars
from path_config import data_path

class sCMD:
	def CMD_gen(self):
		if self.method == 'Vorbin' or self.method == 'kde':
			subprocess.run(['./TestCMDPAR.sh', str(self.mc_num), str(self.mc_num), str(self.pdmf), str(self.binary), '4000000', str(self.age), str(self.GC_name),str(self.feh), '2.0'])
		elif self.method == 'KS2d':
			subprocess.run(['./TestCMDPAR.sh', str(self.mc_num), str(self.mc_num), str(self.pdmf), str(self.binary), '4000000', str(self.age), str(self.GC_name),str(self.feh), '3.0'])

	def CMD_check(self):
		check_time = 0
		exist = False
		while check_time < 3 and exist == False:
			check_time += 1
			if os.path.exists("{}/mc{}.a{}".format(self.outcmd_path,self.mc_num, self.age)):
				exist = True
				break
			else:
				#sleep for 10 second
				time.sleep(10)
		if exist == False:
			self.CMD_gen()
	
	def run_vorbin(self):
		if self.method == 'Vorbin':
			vorbin_chi2.chi2(self.GC_name, self.mc_num, self.age,UniSN=True, write_vorbin=False)
		elif self.method == 'KS2d':
			vorbin_chi2.KS_2d(self.GC_name,self. mc_num, self.age)
		elif self.method == 'kde':
			vorbin_chi2.kde(self.GC_name,self. mc_num, self.age)
		else:
			raise Exception('Cannot find the method to be used to analyse CMD')

	def rmCMD(self):
		#remove searched sCMD file
		file = "mc{}.a{}".format(self.mc_num, self.age)
		path = os.path.join(self.outcmd_path, file)
		if os.path.exists(path) == True:
			os.remove(path)

	def __init__(self, GC_name, mc_num, age, method):
		self.GC_name = str(GC_name)
		#define global variables
		self.outcmd_path = "{}{}/simulateCMD/outcmd".format(data_path, GC_name)
		self.method = str(method)
		#define globar variables used to generate sCMD
		self.feh, self.binary, self.pdmf = sCMD_vars(GC_name)
		self.mc_num = str(mc_num)
		if float(age) < 10000:
			self.age = '0'+str(age)
		else:
			self.age = str(age)
		self.outchi2_path = "{}{}/outchi2/chi2_a{}_mc{}".format(data_path,GC_name,self.age,self.mc_num)
		#if the chi2 file already exist, skip
		if os.path.exists(self.outchi2_path) == False:
			#time it
			start_time = time.time()
			#generate CMD
			self.CMD_gen()
			self.CMD_check()
			self.run_vorbin()
			self.rmCMD()
			print("Done MC{} Age{}".format(self.mc_num, self.age))
			print("--- %s seconds ---" % (time.time() - start_time))
		else:
			self.rmCMD()