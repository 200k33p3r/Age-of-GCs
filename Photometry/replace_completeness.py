#This piece of code is used to replace the completeness with the
#completeness determined from the AS test

import sys
import os
GC_name=str(sys.argv[1])
fidanka_path="{}_fidanka_inputfiles".format(GC_name)
AS_path="{}/{}_inputfiles".format(os.getcwd(),GC_name)
os.chdir(fidanka_path)
for filename in os.listdir():
    if filename[1:4] == 'err':
        if filename[4:] == 'Boundarys.dat':
            target = open(filename, "r+")
            lines = target.readlines()
            AS_file = open("{}/{}".format(AS_path,filename), "r+")
            AS_lines = AS_file.readlines()
            target.close()
            AS_file.close()
            for i in range(3,len(lines)):
                lines[i] = lines[i][:-16] + AS_lines[i][-16:]
            out = open(filename, 'w')
            out.writelines(lines)
            out.close()