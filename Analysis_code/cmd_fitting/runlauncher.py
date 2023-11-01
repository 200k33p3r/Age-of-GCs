#import runsCMD
import sys
import os
import numpy as np
import time

GC_name = str(sys.argv[1])
#mc_num = str(sys.argv[2])
#age = str(sys.argv[3])
#method = str(sys.argv[4])

order_num = str(sys.argv[2])
ages=np.linspace(8000,16000,41)
mc_num = str(20000 + int(order_num)//2)

method = str(sys.argv[3])

#print("GC_name is {}, Age is: {}, MC_num is: {}".format(GC_name, age, mc_num))
print("GC_name is {}, MC_num is: {}".format(GC_name, mc_num))

tmp = "/scratch/" + str(int(time.time())) + str(np.random.randint(1000000000,high=9999999999))
os.mkdir(tmp)
#copy files to tmp 
filelist = ['runsCMD.py', 'vorbin_chi2.py','TestCMDPAR.sh']
for file in filelist:
    os.system("cp {} {}".format(file,tmp))
os.chdir(tmp)
import runsCMD
# runsCMD.sCMD(GC_name,mc_num,age,method)
if int(order_num) % 2 == 0:
    for age in ages[:20]:
        runsCMD.sCMD(GC_name,mc_num,int(age),method)
else:
    for age in ages[20:]:
        runsCMD.sCMD(GC_name,mc_num,int(age),method)
