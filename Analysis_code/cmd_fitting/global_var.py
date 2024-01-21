def sCMD_vars(GC_name):
    if GC_name == 'M55':
        feh=190
        #binary from Milone 2012 A&A 540, A16 (2012)
        binary=0.04
        #binaries are already included in the phonetric error
        #self.binary=0.00
        pdmf = -0.83
    elif GC_name == 'NGC3201':
        feh=148
        binary=0.061
        pdmf = -1.22
    elif GC_name == 'M15':
        feh=227
        binary=0.013
        pdmf=-0.99
    elif GC_name == 'M30':
        #pdmf from Ebrahimi et al (2020)
        feh=210
        binary=0.012
        pdmf=-0.80
    return feh, binary, pdmf

def define_range(GC_name):
    if GC_name == 'M92':
        feh = 230
    elif GC_name == 'M55':
        # dm_max = 14.40
        # dm_min = 13.40
        # red_max = 0.20
        # red_min = 0.00
        dm_max = 14.1
        dm_min = 13.8
        red_max = 0.15
        red_min = 0.08
        feh = 190
    elif GC_name == 'NGC3201':
        dm_max = 14.3
        dm_min = 14.0
        red_max = 0.30
        red_min = 0.15
        feh=148
    elif GC_name == 'M15':
        dm_max = 15.8
        dm_min = 15.3
        red_max = 0.15
        red_min = 0.08
        feh=227
    return feh, dm_max, dm_min, red_max, red_min

def define_N_true_obs(GC_name):
    if GC_name == 'M15':
        N_true_obs = 91944
    return N_true_obs

def define_N_phot(GC_name):
    if GC_name == 'M92':
        N_phot = 80
    elif GC_name == 'NGC3201':
        N_phot = 120
    elif GC_name == 'M55':
        N_phot = 120
    elif GC_name == 'M15':
        N_phot = 120  
    return N_phot