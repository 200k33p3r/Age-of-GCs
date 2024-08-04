import DSEP_code.MC_var_generator_sobol.sobol_mcvar as sobol_mcvar
import os

outpath = os.getcwd()
mc_start, mc_end = 10000, 20000
feh = 0.0
sigma_feh = 0.0
afe = 0.0
sigma_afe = 0.0
sobol_mcvar.genetatemcvar(outpath, mc_start, mc_end, feh = feh, sigma_feh = sigma_feh, afe = afe, sigma_afe = sigma_afe)