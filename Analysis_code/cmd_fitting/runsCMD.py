import vorbin_chi2
import subprocess
import os
import time

class sCMD:
	def CMD_gen(self):
		subprocess.run(['./TestCMDPAR.sh', str(self.mc_num), str(self.mc_num), str(self.pdmf), str(self.binary), '4000000', str(self.age), str(self.GC_name),str(self.feh)])

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
		if self.method == 'Vorbin':
			vorbin_chi2.chi2(self.GC_name, self.mc_num, self.age,UniSN=False, write_vorbin=False)
		elif self.method == 'KS2d':
			vorbin_chi2.KS_2d(self.GC_name,self. mc_num, self.age)
		else:
			raise Exception('Cannot find the method to be used to analyse CMD')

	def rmCMD(self):
		#remove searched sCMD file
		file = "mc{}.a{}".format(self.mc_num, self.age)
		path = os.path.join(self.outcmd_path, file)
		os.remove(path)

	def __init__(self, GC_name, mc_num, age, method):
		#define global variables
		self.outcmd_path = "/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/{}/simulateCMD/outcmd".format(GC_name)
		self.GC_name = str(GC_name)
		self.method = str(method)
		if GC_name == 'M55':
			self.feh=190
			#binary from Milone 2012 A&A 540, A16 (2012)
			self.binary=0.04
			#binaries are already included in the phonetric error
			#self.binary=0.00
			self.pdmf = -0.83
		elif GC_name == 'NGC3201':
			self.feh=148
			self.binary=0.061
			self.pdmf = -1.22
		elif GC_name == 'M15':
			self.feh=227
			self.binary=0.013
			self.pdmf=-0.99
		elif GC_name == 'M30':
			#pdmf from Ebrahimi et al (2020)
			self.feh=210
			self.binary=0.012
			self.pdmf=-0.80
		self.mc_num = str(mc_num)
		if float(age) < 10000:
			self.age = '0'+str(age)
		else:
			self.age = str(age)
		self.outchi2_path = "/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/{}/outchi2/chi2_a{}_mc{}".format(GC_name,self.age,self.mc_num)
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