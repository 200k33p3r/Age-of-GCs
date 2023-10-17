#import runsCMD
import sys
import os
import numpy as np
import time

GC_name = str(sys.argv[1])
mc_num = str(sys.argv[2])
age = str(sys.argv[3])
method = str(sys.argv[4])
print("GC_name is {}, Age is: {}, MC_num is: {}".format(GC_name, age, mc_num))

tmp = "/scratch/08819/mying/" + str(int(time.time())) + str(np.random.randint(1000000000,high=9999999999))
os.mkdir(tmp)
#copy files to tmp 
filelist = ['runsCMD.py', 'vorbin_chi2.py','TestCMDPAR.sh']
for file in filelist:
    os.system("cp {} {}".format(file,tmp))
os.chdir(tmp)
import runsCMD
runsCMD.sCMD(GC_name,mc_num,age,method)
