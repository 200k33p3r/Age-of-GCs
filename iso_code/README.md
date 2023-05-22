# New Isochrone Construction Code

I've been looking at creating synthetic colour-magnitude diagrams (CMD), and realized that we have to modify how we are creating isochrones from the DSEP stellar evolution tracks. We need a much finer mass spacing in the isochrones. Isochrone generation for synthetic CMD production is now a 3 step process, as outlined below. I've created a sample csh script which does all three steps, called `isomc.sh`.    You will need to modify `isomc.sh` so that it know the location of the various input/output directories. As with the previous version of the code, the subdirectores `/phx` and `/vdb` need to exist in the run directory, as these subdirectories contain the bolometric corrections files.  Otherwise, I think it will work fine with the your Monte Carlo DSEP track files, as I tested it on a bunch of your files which I downloaded from discovery. 

## Step 1

Compile `isomc2.f` using: `gfortran  -std=legacy isomc2.f -o isomc2` . 

Run `isomc2` , which I believe is identical to the code you already have.  However, we have to modify the input namelist file to tell the code to generate a much finer mass grid.  The following lines need to be added to the input namelist (called `isomc.nml`  in my sample csh script).   Note that your namelist file has many more ages, so likely best if you just add the lines below to your existing namelist file. 

```
  QT(1) = 300
  QT(2) = 100
  QT(3) = 80
  QT(4) = 30
  QT(5) = 50
  QT(6) = 35
  QT(7) = 10
  QT(8) = 8
  QT(9) = 8  
  QT(10) = 8
  QT(11) = 8
  QT(12) = 0
```

## Step 2
Compile `lowmass2.f`  using :  `gfortran lowmass2.f -o lowmass2` 

Run `lowmass2` which is a modified version of the lowmass isochrone code you have been running.  This code replaces the low mass points in the input isochorne (created with `isomc2` with points derived from the low mass DSEP runs which included used FreeEOS.  In addition to the input isochrone and track files (set in `isomc.sh` ) the code requires that you enter at the command line the number of low mass tracks, the number of ages in the isochrone file and the minimum mass which is NOT being replaced in the isochorne file.  The first two numbers are already in the input namelist file to isomc2, so see line 73 in `isomc.sh`  for how to set these numbers automatically.

## Step 3
Compile `prepiso.f` using : `gfortran prepiso.f -o prepiso`

Run `prepiso` which reads in the output from `lowmass2` and outputs a file for use by the synthetic CMD program.  This includes doing a smooth interpoaltion among the existing masss points, and outputing a file which just contains mass, and magnitudes/colours in F606W and F814W.

