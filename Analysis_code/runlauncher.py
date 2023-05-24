import runsCMD
import sys

GC_name = str(sys.argv[1])
mc_num = str(sys.argv[2])
age = str(sys.argv[3])
print("GC_name is {}, Age is: {}, MC_num is: {}".format(GC_name, age, mc_num))
runsCMD.sCMD(GC_name,mc_num,age)
