# -*- coding: utf-8 -*-
"""
User settings parameters

@author: Prasanna Sritharan
"""



'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''


'''
UserSettings:
    Parent user settings class for C3D processing pipeline.
'''
class UserSettings():
    def __init__(self):
        
        
        # ******************************
        # GENERAL SETTINGS        
        
        # data folders
        self.rootpath = r"C:\Users\Owner\Documents\data\FORCe"
        self.infolder = [r"inputdatabase\ForceMaster_LTUControls\3.2. OpenSim", r"C:\Users\Owner\Documents\data\FORCe\inputdatabase\ForceMaster_LTUSymptomatics\3.2. OpenSim_C3D_Good"]
        self.outfolder = "outputdatabase"
        self.trialgroupfolders = [""]    


        # ******************************
        # C3D DATA PROCESSING  
        
        # Add generic parameters here
        
        
        # ******************************
        # OPENSIM PARAMETERS
        
        # OpenSim log file
        self.logfilepath = r"C:\Users\Owner\Documents\projects\force-moco\python"
        self.triallogfolder = "log"
        self.logfile = "opensim.log"
        
        # OpenSim reference model
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\force-moco\python\opensim-models"
        self.refmodelfile = "LASEM_FORCE_ReferenceModel_PelvisUnclamped.osim"

        # OpenSim setup files and folders
        self.refsetuppath = r"C:\Users\Owner\Documents\projects\force-moco\python\opensim-setup"
        self.refsetupscale = "LASEM_FORCE_Setup_Scale.xml"
        self.refsetupik = "LASEM_FORCE_Setup_IK.xml"
        self.refsetupid = "LASEM_FORCE_Setup_ID.xml"
        self.refsetupso = "LASEM_FORCE_Setup_Analysis.xml"
        self.refsetuprra = "LASEM_FORCE_Setup_RRA.xml"
        self.refsetupcmc = "LASEM_FORCE_Setup_CMC.xml"   
        
        # OpenSim analysis codes
        self.scalecode = "scale"
        self.ikcode = "ik"
        self.idcode = "id"
        self.socode = "so"
        self.rracode = "rra"
        self.cmccode = "cmc"
        self.jrcode = "jr"
        
        
        
'''
FORCESettings_SDP(UserSettings):
    UserSettings for LASEM FORCE Project: SDP
'''
class FORCESettings_SDP(UserSettings):
    def __init__(self):
        
        
        # inherit parent attributes
        super(FORCESettings_SDP, self).__init__()
        
                        
        # ******************************
        # GENERAL SETTINGS
        
        # project
        self.project = "FORCE_SDP"

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "force_sdp_results_all_"
        
        # meta data file
        self.metadatafile = self.project + ".pkl"        
        
        # file prefixes
        self.subjprefix = "FAILT"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static01"
        self.staticfpchannel = "Fz1"
        
        # file suffixes based on task
        self.trialprefixes = {}
        self.trialprefixes["sdp"] = ["SDP"]
                        
        # file name format regex pattern:
        #   (subjprefix)(num code)_(trialprefix)(num code)
        self.fnpat = "(FAILTCRT|FAILT)(\d+)_(STATIC|SDP)(\d*)([A-Za-z]*)" #"(FAILTCRT|FAILT)(\d+)_([A-Za-z]+)(\d+)([A-Za-z]*)"
        self.tasktoknum = 3   # the token + 1 that represents the task name/type
        
        # output samples
        self.samples = 101
       
        
        # ******************************
        # C3D DATA PROCESSING     
       
        # marker data filter (set cutoff to -1 if not required)
        self.marker_filter_butter_order = 4
        self.marker_filter_cutoff = -1
       
        # force plate data filter (set cutoff to -1 if not required)
        self.fp_filter_butter_order = 4
        self.fp_filter_cutoff = 15
        self.fp_filter_threshold = 15
        self.fp_smooth_transitions = False
        self.fp_smooth_cop_fixed_offset = 25   # required but not currently used
        self.fp_smooth_window = 20
        
        
        # ******************************
        # OPENSIM PARAMETERS
        
        # OpenSim additional files
        self.additionalfilesfolder = "SDP"
        self.refexternalloads = "LASEM_FORCE_ExternalLoads.xml"
        self.refreserveactuators = "LASEM_FORCE_Reserve_Actuators.xml"
        self.refrraactuators = "LASEM_FORCE_RRA_Actuators_RUN.xml"
        self.refrratasks = "LASEM_FORCE_RRA_Tasks_RUN.xml"
        self.refcmcactuators = "LASEM_FORCEL_CMC_Actuators_RUN.xml"
        self.refcmctasks = "LASEM_FORCE_CMC_Tasks_RUN.xml"
        self.refcmccontrolconstraints = "LASEM_FORCE_CMC_ControlConstraints_RUN.xml"
        
        # OpenSim Scale parameters
        #self.fom_scalefactor = {}
        #self.fom_scalefactor["all"] = 3.0
        #self.lom_lmt_scalefactor = {}
        #self.lom_lmt_scalefactor["all"] = 1.1
        
        # OpenSim IK parameters
        self.kinematics_filter_butter_order = 4
        self.kinematics_filter_cutoff = 6
        
        # OpenSim RRA parameters
        self.update_mass = True
        self.rraiter = 2   
        self.rra_start_time_offset = -0.03  # to enable CMC initalisation
        self.rra_end_time_offset = 0.03     # slightly wider than CMC end time       
        self.prescribe_upper_body_motion = True
        self.prescribe_coord_list = ["lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex_r", "arm_add_r", "arm_rot_r", "elbow_flex_r", "pro_sup_r", "wrist_flex_r", "wrist_dev_r", "arm_flex_l", "arm_add_l", "arm_rot_l", "elbow_flex_l", "pro_sup_l", "wrist_flex_l", "wrist_dev_l"]

        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_rra_kinematics = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.03  # to enable CMC initalisation
        self.cmc_end_time_offset = 0.0
        
        # OpenSim JR parameters
        self.jr_joints = {}
        self.jr_joints["all"] = ["child", "child"]
        self.jr_use_cmc_forces = False


        # ******************************
        # OPENSIM RESULTS
        
        # left leg flip columns (incl. time)
        self.results_flip = {}
        self.results_flip["ik"] = [2, 3, 6, 24, 25]
        self.results_flip["id"] = [2, 3, 6, 14, 15]
        self.results_flip["so"] = []
        self.results_flip["rra"] = []
        self.results_flip["cmc"] = []
        self.results_flip["jr"] = []  
        
        # foot columns (incl. time): R, L
        self.results_columns = {}
        self.results_columns["ik"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], 
                                      [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 33, 34, 35, 36, 37, 38, 39]]
        self.results_columns["id"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 20, 21, 22, 26, 28, 30, 32, 34, 36, 37], 
                                      [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 18, 19, 23, 24, 25, 27, 29, 31, 33, 35, 38, 39]]
        self.results_columns["so"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90],
                                      [0, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 91, 92, 93, 94, 95, 96, 97]]
        self.results_columns["rra"] = []
        self.results_columns["cmc"] = []
        self.results_columns["jr"] = []   
        
        # headers
        self.results_headers = {}
        self.results_headers["ik"] = ["time", "pelvis_tilt", "pelvis_list", "pelvis_rotation", "pelvis_tx", "pelvis_ty", "pelvis_tz", "hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "knee_angle_beta", "ankle_angle", "subtalar_angle", "mtp_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex", "arm_add", "arm_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["id"] = ["time", "pelvis_tilt_moment", "pelvis_list_moment", "pelvis_rotation_moment", "pelvis_tx_force", "pelvis_ty_force", "pelvis_tz_force", "hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment", "knee_angle_moment", "knee_angle_beta_force", "arm_flex_moment", "arm_add_moment", "arm_rot_moment", "ankle_angle_moment", "elbow_flex_moment", "subtalar_angle_moment", "pro_sup_moment", "mtp_angle_moment", "wrist_flex_moment", "wrist_dev_moment"]
        self.results_headers["so"] = ["time", "addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx", "bflh", "bfsh", "edl", "ehl", "fdl", "fhl", "gaslat", "gasmed", "glmax1", "glmax2", "glmax3", "glmed1", "glmed2", "glmed3", "glmin1", "glmin2", "glmin3", "grac", "iliacus", "perbrev", "perlong", "piri", "psoas", "recfem", "sart", "semimem", "semiten", "soleus", "tfl", "tibant", "tibpost", "vasint", "vaslat", "vasmed", "lumbar_ext", "lumbar_bend", "lumbar_rot", "shoulder_flex", "shoulder_add", "shoulder_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["rra"] = []
        self.results_headers["cmc"] = []
        self.results_headers["jr"] = []          
        

        
'''
TRAILSettings_RUN(UserSettings):
    UserSettings for LASEM TRAIL Project: RUN
'''
class TRAILSettings_RUN(UserSettings):
    def __init__(self):
        
        
        # inherit parent attributes
        super(TRAILSettings_RUN, self).__init__()
        
                        
        # ******************************
        # GENERAL SETTINGS
        
        # project
        self.project = "TRAIL_RUN"

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "trail_run_opensim_results_all_"
        
        # meta data file
        self.metadatafile = self.project + ".pkl"        
        
        # file prefixes
        self.subjprefix = "TRAIL_"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static_01"
        self.staticfpchannel = "Force.Fz3"
        
        # file suffixes based on task
        self.trialprefixes = {}
        self.trialprefixes["run_stance"] = ["EP", "FAST"]
        self.trialprefixes["run_stridecycle"] = ["EP", "FAST"]
        self.trialprefixes["run_stance_ep"] = ["EP"]
        self.trialprefixes["run_stridecycle_ep"] = ["EP"]
        self.trialprefixes["run_stance_fast"] = ["FAST"]
        self.trialprefixes["run_stridecycle_fast"] = ["FAST"]
                        
        # file name format regex pattern:
        #   (subjprefix)_(num code)_(trialprefix)_(alphanum code)
        self.fnpat = "TRAIL_\d+_(\w+)_\w+"
        self.tasktoknum = 1   # the token + 1 that represents the task name/type
        
        # output samples
        self.samples = 101
       
        
        # ******************************
        # C3D DATA PROCESSING     
       
        # marker data filter (set cutoff to -1 if not required)
        self.marker_filter_butter_order = 4
        self.marker_filter_cutoff = -1
       
        # force plate data filter (set cutoff to -1 if not required)
        self.fp_filter_butter_order = 4
        self.fp_filter_cutoff = 15
        self.fp_filter_threshold = 15
        self.fp_smooth_transitions = True
        self.fp_smooth_cop_fixed_offset = 25   # required but not currently used
        self.fp_smooth_window = 20
        
        
        # ******************************
        # OPENSIM PARAMETERS
        
        # OpenSim additional files
        self.additionalfilesfolder = "RUN"
        self.refexternalloads = "LASEM_TRAIL_ExternalLoads.xml"
        self.refreserveactuators = "LASEM_TRAIL_Reserve_Actuators_WithUpper.xml"
        self.refrraactuators = "LASEM_TRAIL_RRA_Actuators_RUN.xml"
        self.refrratasks = "LASEM_TRAIL_RRA_Tasks_RUN.xml"
        self.refcmcactuators = "LASEM_TRAIL_CMC_Actuators_RUN.xml"
        self.refcmctasks = "LASEM_TRAIL_CMC_Tasks_RUN.xml"
        self.refcmccontrolconstraints = "LASEM_TRAIL_CMC_ControlConstraints_RUN.xml"
        
        # OpenSim Scale parameters
        self.fom_scalefactor = {}
        self.fom_scalefactor["all"] = 3.0
        #self.lom_lmt_scalefactor = {}
        #self.lom_lmt_scalefactor["all"] = 1.1
        
        # OpenSim IK parameters
        self.kinematics_filter_butter_order = 4
        self.kinematics_filter_cutoff = 6
        
        # OpenSim RRA parameters
        self.update_mass = True
        self.rraiter = 2   
        self.rra_start_time_offset = -0.03  # to enable CMC initalisation
        self.rra_end_time_offset = 0.03     # slightly wider than CMC end time       
        self.prescribe_upper_body_motion = True
        self.prescribe_coord_list = ["lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex_r", "arm_add_r", "arm_rot_r", "elbow_flex_r", "pro_sup_r", "wrist_flex_r", "wrist_dev_r", "arm_flex_l", "arm_add_l", "arm_rot_l", "elbow_flex_l", "pro_sup_l", "wrist_flex_l", "wrist_dev_l"]

        # OpenSim CMC parameters
        self.use_rra_model = True
        self.use_rra_kinematics = True
        self.use_fast_target = True
        self.cmc_start_time_offset = -0.03  # to enable CMC initalisation
        self.cmc_end_time_offset = 0.0
        
        # OpenSim JR parameters
        self.jr_joints = {}
        self.jr_joints["all"] = ["child", "child"]
        self.jr_use_cmc_forces = False        
        
        
        # ******************************
        # OPENSIM RESULTS
        
        # left leg flip columns (incl. time)
        self.results_flip = {}
        self.results_flip["ik"] = [3, 4, 7, 25, 26]
        self.results_flip["id"] = [3, 4, 7, 15, 16]
        self.results_flip["so"] = []
        self.results_flip["rra"] = []
        self.results_flip["cmc"] = []
        self.results_flip["jr"] = []  
                
        # foot columns (incl. time): R, L
        self.results_columns = {}
        self.results_columns["ik"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], 
                                      [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 33, 34, 35, 36, 37, 38, 39]]
        self.results_columns["id"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 20, 21, 22, 26, 28, 30, 32, 34, 36, 37], 
                                      [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 18, 19, 23, 24, 25, 27, 29, 31, 33, 35, 38, 39]]
        self.results_columns["so"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90],
                                      [0, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 91, 92, 93, 94, 95, 96, 97]]
        self.results_columns["rra"] = []
        self.results_columns["cmc"] = []
        self.results_columns["jr"] = []   
        
        # headers
        self.results_headers = {}
        self.results_headers["ik"] = ["time", "pelvis_tilt", "pelvis_list", "pelvis_rotation", "pelvis_tx", "pelvis_ty", "pelvis_tz", "hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "knee_angle_beta", "ankle_angle", "subtalar_angle", "mtp_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex", "arm_add", "arm_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["id"] = ["time", "pelvis_tilt_moment", "pelvis_list_moment", "pelvis_rotation_moment", "pelvis_tx_force", "pelvis_ty_force", "pelvis_tz_force", "hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment", "knee_angle_moment", "knee_angle_beta_force", "arm_flex_moment", "arm_add_moment", "arm_rot_moment", "ankle_angle_moment", "elbow_flex_moment", "subtalar_angle_moment", "pro_sup_moment", "mtp_angle_moment", "wrist_flex_moment", "wrist_dev_moment"]
        self.results_headers["so"] = ["time", "addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx", "bflh", "bfsh", "edl", "ehl", "fdl", "fhl", "gaslat", "gasmed", "glmax1", "glmax2", "glmax3", "glmed1", "glmed2", "glmed3", "glmin1", "glmin2", "glmin3", "grac", "iliacus", "perbrev", "perlong", "piri", "psoas", "recfem", "sart", "semimem", "semiten", "soleus", "tfl", "tibant", "tibpost", "vasint", "vaslat", "vasmed", "lumbar_ext", "lumbar_bend", "lumbar_rot", "shoulder_flex", "shoulder_add", "shoulder_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["rra"] = []
        self.results_headers["cmc"] = []
        self.results_headers["jr"] = []   