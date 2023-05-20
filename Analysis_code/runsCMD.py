import vorbin_chi2
import subprocess
import os
import time

class sCMD:
	def CMD_gen(self):
		subprocess.run(['./TestCMDPAR.sh', str(self.mc_num), str(self.mc_num), '-1.02', '0.02', '4000000', str(self.age), str(self.GC_name)])

	def CMD_check(self):
		check_time = 0
		exist = False
		while check_time < 10 and exist == False:
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
		vorbin_chi2.chi2(self.GC_name, self.mc_num, self.age)

	def rmCMD(self):
		#remove searched sCMD file
		file = "mc{}.a{}".format(self.mc_num, self.age)
		path = os.path.join(self.outcmd_path, file)
		os.remove(path)

	def __init__(self, GC_name, mc_num, age):
		#define global variables
		self.outcmd_path = "/work2/08819/mying/{}/simulateCMD/outcmd".format(GC_name)
		self.GC_name = str(GC_name)
		self.mc_num = str(mc_num)
		self.age = str(age)
		#time it
		start_time = time.time()
		#generate CMD
		self.CMD_gen()
		self.CMD_check()
		self.run_vorbin()
		self.rmCMD()
		print("Done MC{} Age{}".format(self.mc_num, self.age))
		print("--- %s seconds ---" % (time.time() - start_time))