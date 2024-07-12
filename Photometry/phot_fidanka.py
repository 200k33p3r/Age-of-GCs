#This code is designed to determine the photometric error for HST ACS globular cluster data
#in case the photometric uncertainty estimated from the artifical star test is not sufficient
#to explain the width in the CMD.
#This is likely to happen in GCs which does not suffer from crowding.
#The result file is an sudo artificial star test result that can be feeded into the photometry code

import numpy as np
import pandas as pd
import os
import shapely.geometry as geom
import sys

class phot:
    #read in fitstars and fiducial isochrone generated from the fidanka package
    def read_input(self, data_path, iso_path):
        fitstars = pd.read_csv(data_path)
        iso = pd.read_csv(iso_path)
        v_star = fitstars['v'].values
        i_star = fitstars['i'].values
        v_iso = iso['v'].values
        i_iso = iso['v'].values - iso['vi'].values
        x_star = fitstars['x'].values
        y_star = fitstars['y'].values
        return v_star, i_star, v_iso, i_iso, x_star, y_star
    
    #Even though the artificial star test is not capable of estimating the photometric
    #error in certain cases, it is still helpfuly in determine the completeness.
    #In this case, we adopt the completeness determined from AS test
    # def read_completeness(self, input_file_path):
    #     #do v first
    #     Verr = pd.read_csv("{}/VerrBoundart.dat",skiprows=3, names=['outboundRad','minMagBound','Nstar','Completenss'])
    #     #do i
    #     Ierr = pd.read_csv("{}/IerrBoundart.dat",skiprows=3, names=['outboundRad','minMagBound','Nstar','Completenss'])

    
    def determine_error(self, v_star, i_star, v_iso, i_iso):
        coords = np.vstack((v_iso,i_iso)).T
        line = geom.LineString(coords)
        out_v = np.zeros(len(v_star))
        out_i = np.zeros(len(i_star))
        for i in range(len(v_star)):
            point = geom.Point(np.array([v_star[i], i_star[i]]))
            point_on_line = line.interpolate(line.project(point))
            out_v[i] = point_on_line.x
            out_i[i] = point_on_line.y
        return out_v, out_i
    
    #write the fiducial AS test
    def write_out(self, fiducial_AS_path, v_star, i_star, out_v, out_i, x_star, y_star):
        data = np.vstack((x_star, y_star, v_star - self.V_diff, i_star - self.I_diff, x_star, y_star, out_v - self.V_diff, out_i - self.I_diff)).T
        df = pd.DataFrame(data=data, columns=['InputX','InputY','InputF606W','InputF814W','OutputX','OutputY','OutputF606W','OutputF814W'])
        df.to_csv(fiducial_AS_path,index=False,sep=' ')
    
    def __init__(self,GC_name):
        if GC_name == 'M55':
            self.V_diff = 30.98
            self.I_diff = self.V_diff-0.737
        elif GC_name == 'NGC3201':
            self.V_diff = 31.4
            self.I_diff = self.V_diff-0.897
        elif GC_name == 'M15':
            self.V_diff = 31.638
            self.I_diff = self.V_diff-0.728
        elif GC_name == 'M30':
            self.V_diff = 31.70
            self.I_diff = self.V_diff-0.885
        #define all the path for read and write
        cwd = os.getcwd()
        repo_path = os.path.abspath(os.path.join(cwd, os.pardir))
        data_path = "{}/{}_data".format(repo_path, GC_name)
        obs_data_path = "{}/{}_fitstars.dat".format(data_path,GC_name)
        iso_path = "{}/fiducial_lines.csv".format(data_path)
        fiducial_AS_path = "{}/{}_fiducial_artstars.dat".format(data_path,GC_name)

        #read files
        v_star, i_star, v_iso, i_iso, x_star, y_star = self.read_input(obs_data_path, iso_path)
        #determine error
        out_v, out_i = self.determine_error(v_star, i_star, v_iso, i_iso)
        #write output
        self.write_out(fiducial_AS_path, v_star, i_star, out_v, out_i, x_star, y_star)

#read input and run code
GC_name = str(sys.argv[1])
phot(GC_name)