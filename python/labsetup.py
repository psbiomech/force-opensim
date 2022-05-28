# -*- coding: utf-8 -*-
"""
Lab set ups for C3D extract

@author: Prasanna Sritharan
"""


import numpy as np




'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''


'''
LabKey():
    Template class for lab set up info.
'''
class LabKey():
    def __init__(self):
        
        # lab name
        self.lab_name = "NoName"
        
        # force plates
        self.fp_type = -1
        self.n_fp = -1
        self.n_fp_channels = -1
        self.fp_used = 0
        self.fp_dict_name_prefix = ""
        self.fp_channel_prefixes = []
        self.fp_channel_suffixes = []
        
        # transform vectors
        self.transform_lab_to_opensim_xdir = []
        self.transform_fp_to_lab = []
        
        # transform matrices from transform vectors
        self.transform_mat_lab_to_opensim = [np.zeros([3,3])]*6
        self.transform_mat_fp_to_lab = np.zeros([3,3])
        
        # markers
        self.marker_list = []
        self.offset_marker = []
        
        return None



'''
labKeyLasemForceSDP(LabKey):
    LabKey for LASEM FORCe step-down-and-pivot (SDP).
'''
class LabKeyLasemForceSDP(LabKey):
    def __init__(self):
        
        self.lab_name = "lasem_trail"    
    
        # force plates
        self.fp_type = 2
        self.n_fp = 4
        self.n_fp_channels = 6
        self.fp_used = [1, 2, 3]  # plates used for dynamic trials
        self.fp_used_static = [1, 2]  # plates used for static trials
        self.fp_dict_name_prefix = "FP"
        self.fp_channel_prefixes = ["Fx", "Fy", "Fz", "Mx", "My", "Mz"]
        self.fp_channel_suffixes = ["", "", "", "", "", ""]
        
        # transform vectors
        self.transform_lab_to_opensim_xdir = [1, 3, -2]
        self.transform_fp_to_lab = [-1, 2, -3]  # temporary, better to take this from C3D file as force plates can be moved around
        
        # transform matrices from transform vectors
        self.transform_mat_lab_to_opensim = create_transform_set_lab_to_opensim(self.transform_lab_to_opensim_xdir)
        self.transform_mat_fp_to_lab = create_transform_matrix(self.transform_fp_to_lab)
                
        # markers
        self.marker_list = ["RSH", "LSH", "C7", "LELB", "LWR", "RELB", "RWR", "RASI", "LASI", "SACR", "P1", "P2", "LTHLP", "LTHLD", "LTHAP", "LTHAD", "LLEPI", "LPAT", "LTIAP", "LTIAD", "LTILAT", "LLMAL", "LHEEL", "LMFS", "LMFL", "LP5MT", "LP1MT", "LTOE", "RTHLP", "RTHLD", "RTHAP", "RTHAD", "RLEPI", "RPAT", "RTIAP", "RTIAD", "RTILAT", "RLMAL", "RHEEL", "RMFS", "RMFL", "RP5MT", "RP1MT", "RTOE"]
        self.offset_marker = "SACR"

        return None



'''
labKeyLasemTrail(LabKey):
    LabKey for LASEM TRAIL project.
'''
class LabKeyLasemTrail(LabKey):
    def __init__(self):
        
        self.lab_name = "lasem_trail"    
    
        # force plates
        self.fp_type = 2
        self.n_fp = 4
        self.n_fp_channels = 6
        self.fp_used = [3, 4]  # plates used for dynamic trials
        self.fp_used_static = [3, 4]  # plates used for static trials
        self.fp_dict_name_prefix = "FP"
        self.fp_channel_prefixes = ["Force.Fx", "Force.Fy", "Force.Fz", "Moment.Mx", "Moment.My", "Moment.Mz"]
        self.fp_channel_suffixes = ["", "", "", "", "", ""]
        
        # transform vectors
        self.transform_lab_to_opensim_xdir = [1, 3, -2]
        self.transform_fp_to_lab = [-1, 2, -3]  # temporary, better to take this from C3D file as force plates can be moved around
        
        # transform matrices from transform vectors
        self.transform_mat_lab_to_opensim = create_transform_set_lab_to_opensim(self.transform_lab_to_opensim_xdir)
        self.transform_mat_fp_to_lab = create_transform_matrix(self.transform_fp_to_lab)
                
        # markers
        self.marker_list = ["RSH", "LSH", "C7", "LELB", "LWR", "RELB", "RWR", "RASI", "LASI", "SACR", "P1", "P2", "LTHLP", "LTHLD", "LTHAP", "LTHAD", "LLEPI", "LPAT", "LTIAP", "LTIAD", "LTILAT", "LLMAL", "LHEEL", "LMFS", "LMFL", "LP5MT", "LP1MT", "LTOE", "RTHLP", "RTHLD", "RTHAP", "RTHAD", "RLEPI", "RPAT", "RTIAP", "RTIAD", "RTILAT", "RLMAL", "RHEEL", "RMFS", "RMFL", "RP5MT", "RP1MT", "RTOE"]
        self.offset_marker = "SACR"

        return None



'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
create_transform_set_lab_to_opensim(tvec_x):
    Create lab to opensim rotation matrices from transform vectors. Note: Z 
    transforms are assumed and may need to be respecified depending on the lab 
    set-up as there are 2 possible Z-direction configurations.
'''
def create_transform_set_lab_to_opensim(tvec_x):
    
    # initiliase set of transform matrices for directions: X -X Y -Y Z -Z
    transform_set = [None]*6
    
    # Vicon X ---> OpenSim X
    transform_set[0] = create_transform_matrix(tvec_x)
    
    # Vicon -X ---> OpenSim X
    tvec_nx = np.multiply(tvec_x,[-1, 1, -1])
    transform_set[1] = create_transform_matrix(tvec_nx)
    
    # Vicon Y ---> Opensim X
    tvec_y = [-tvec_x[2], tvec_x[1], tvec_x[0]]
    transform_set[2] = create_transform_matrix(tvec_y)
    
    # Vicon -Y ---> OpenSim X
    tvec_ny = np.multiply(tvec_y,[-1, 1, -1])
    transform_set[3] = create_transform_matrix(tvec_ny)    
    
    # Vicon Z ---> OpenSim X (assumed, change if necessary)
    tvec_z = [tvec_x[1], -tvec_x[2], -tvec_x[0]]
    transform_set[4] = create_transform_matrix(tvec_z)  

    # Vicon -Z ---> OpenSim X (assumed, change if necessary)
    tvec_nz = np.multiply(tvec_z,[-1, 1, -1])
    transform_set[5] = create_transform_matrix(tvec_nz)
    
    return transform_set


    
'''
create_transform_matrix(tvec):
    Create transform matrix from transform vectors
'''
def create_transform_matrix(tvec):
    tmat = np.zeros((3,3))
    for r in range(3):
        tmat[r,abs(tvec[r])-1] = 1*np.sign(tvec[r])   
    return tmat

