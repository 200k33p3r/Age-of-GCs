import numpy as np
mc_num_start = 10000
mc_num_end = 11000
GC_name = 'M55'
ages=np.linspace(8000,16000,41)
param_list = []
for mc_num in range(mc_num_start, mc_num_end):
    for age in ages:
    	param_list.append(str(GC_name), (str(mc_num),str(int(age))))
f = open("inputparams", "w")
for i in range(len(param_list)):
    f.write("python runlauncher.py {} {}\n".format(param_list[i][0], param_list[i][1], param_list[i][2]))
f.close()

