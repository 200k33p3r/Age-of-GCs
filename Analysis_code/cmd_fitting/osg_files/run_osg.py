import vorbin_chi2
import subprocess
import os
import time
#modify the runsCMD.py code to run on OSG

def CMD_gen(mc_num, binary, age, GC_name, feh, pdmf):
    subprocess.run(['./TestCMDPAR.sh', str(mc_num), str(mc_num), str(pdmf), str(binary), '4000000', str(age), str(GC_name),str(feh)])

def CMD_check(outcmd_path,mc_num, age):
    check_time = 0
    exist = False
    while check_time < 10 and exist == False:
        check_time += 1
        if os.path.exists("{}/mc{}.a{}".format(outcmd_path,mc_num, age)):
            exist = True
            break
        else:
            #sleep for 10 second
            time.sleep(10)
    if exist == False:
        CMD_gen()

def run_vorbin(GC_name, mc_num, age):
    vorbin_chi2.chi2(GC_name, mc_num, age)

def rmCMD(outcmd_path,mc_num, age):
    #remove searched sCMD file
    file = "mc{}.a{}".format(mc_num, age)
    path = os.path.join(outcmd_path, file)
    os.remove(path)

def run(GC_name, mc_num, age):
    #define global variables
    outcmd_path = os.getcwd()
    GC_name = str(GC_name)
    if GC_name == 'M55':
        feh=190
        #binary from Milone 2012 A&A 540, A16 (2012)
        #binary=0.04
        #binaries are already included in the phonetric error
        binary=0.00
        #pdmf from Ebrahimi et al 2020 MNRAS
        pdmf = -0.83
    elif GC_name == 'M92':
        binary=0.02
        pdmd=-1.02
    mc_num = str(mc_num)
    if float(age) < 10000:
        age = '0'+str(age)
    else:
        age = str(age)
    #time it
    start_time = time.time()
    #generate CMD
    CMD_gen(mc_num, binary, age, GC_name, feh, pdmf)
    CMD_check(outcmd_path,mc_num, age)
    run_vorbin(GC_name, mc_num, age)
    rmCMD(outcmd_path,mc_num, age)
    print("Done MC{} Age{}".format(mc_num, age))
    print("--- %s seconds ---" % (time.time() - start_time))