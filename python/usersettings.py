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
        self.trialgroupfolders = [""]    


        # ******************************
        # C3D DATA PROCESSING  
        
        # Add generic parameters here
        
        
        # ******************************
        # OPENSIM PARAMETERS
        
        # OpenSim log file
        self.logfilepath = r"C:\Users\Owner\Documents\projects\force-opensim\python"
        self.triallogfolder = "log"
        self.logfile = "opensim.log"
        
        # OpenSim reference model
        self.refmodelpath = r"C:\Users\Owner\Documents\projects\force-opensim\python\opensim-models"
        self.refmodelfile = "LASEM_FORCE_ReferenceModel_PelvisUnclamped.osim"

        # OpenSim setup files and folders
        self.refsetuppath = r"C:\Users\Owner\Documents\projects\force-opensim\python\opensim-setup"
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
        self.bkcode = "bk"
        
        # limb code
        self.leg = ["r", "l"]
        
        # gravity
        self.gravity = 9.81
        


'''
FORCESettings_SLDJ(UserSettings):
    UserSettings for LASEM FORCE Project: SLDJ
'''
class FORCESettings_SLDJ(UserSettings):
    def __init__(self):
        
        
        # inherit parent attributes
        super(FORCESettings_SLDJ, self).__init__()
        
                        
        # ******************************
        # GENERAL SETTINGS
        
        # Output database
        self.outfolder = "outputdatabase_sldj"
        
        # project
        self.project = "FORCE_SLDJ"

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "force_sldj_results_all_trials"
        self.csvdescfileprefix = "force_sldj_results_subject_descriptives"
        
        # meta data file
        self.metadatafile = self.project + ".pkl"        
        
        # file prefixes
        self.subjprefix = "FAILT"
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static02"
        self.staticfpchannel = "Fz1"
        
        # file suffixes based on task
        self.trialprefixes = {}
        self.trialprefixes["sldj"] = ["SLDJ"]
                        
        # file name format regex pattern:
        #   (subjprefix)(num code)_(trialprefix)(num code)
        self.fnpat = "(FAILTCRT|FAILT)(\d+)_(STATIC|SLDJ)(\d*)" #"(FAILTCRT|FAILT)(\d+)_([A-Za-z]+)(\d+)([A-Za-z]*)"
        self.tasktoknum = 3   # the token + 1 that represents the task name/type
        
        # output samples
        self.samples = 101

        # additional participant info
        self.additional_participant_info_file = r"outputdatabase\csvfolder\FORCE-ParticipantData-All.csv"
        
        
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
        
        # OpenSim BK parameters
        self.bk_bodies = ["all"]
        self.bk_output_com = True
        self.bk_output_in_local_frame = False
        self.bk_use_cmc_results = False


        # ******************************
        # OPENSIM RESULTS
        
        # left leg flip columns (index incl. time)
        self.results_flip = {}
        self.results_flip["ik"] = [2, 3, 6, 24, 25]
        self.results_flip["id"] = [2, 3, 6, 14, 15]
        self.results_flip["so"] = []
        self.results_flip["rra"] = []
        self.results_flip["cmc"] = []
        self.results_flip["jr"] = [66, 67, 68, 75, 76, 77, 84, 85, 86, 93, 94, 95, 102, 103, 104, 111, 112, 113, 120, 121, 122]  
        self.results_flip["bk"] = range(1, 136)
        
        # foot columns (index incl. time): R, L
        self.results_columns = {}
        self.results_columns["ik"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], 
                                      [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 33, 34, 35, 36, 37, 38, 39]]
        self.results_columns["id"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 20, 21, 22, 26, 28, 30, 32, 34, 36, 37], 
                                      [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 18, 19, 23, 24, 25, 27, 29, 31, 33, 35, 38, 39]]
        self.results_columns["so"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90],
                                      [0, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 91, 92, 93, 94, 95, 96, 97]]
        self.results_columns["rra"] = []
        self.results_columns["cmc"] = []
        self.results_columns["jr"] = [[0, 10, 11, 12, 13, 14, 15, 19, 20, 21, 22, 23, 24, 28, 29, 30, 31, 32, 33, 37, 38, 39, 40, 41, 42, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 118, 119, 120, 121, 122, 123],
                                      [0, 64, 65, 66, 69, 68, 69, 73, 74, 75, 76, 77, 78, 82, 83, 84, 85, 86, 87, 91, 92, 93, 94, 95, 96, 100, 101, 102, 103, 104, 105, 109, 110, 111, 112, 113, 114, 118, 119, 120, 121, 122, 123]]
        self.results_columns["bk"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 133, 134, 135], 
                                      [0, 1, 2, 3, 4, 5, 6, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135]]
        
        # headers
        self.results_headers = {}
        self.results_headers["ik"] = ["time", "pelvis_tilt", "pelvis_list", "pelvis_rotation", "pelvis_tx", "pelvis_ty", "pelvis_tz", "hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "knee_angle_beta", "ankle_angle", "subtalar_angle", "mtp_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex", "arm_add", "arm_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["id"] = ["time", "pelvis_tilt_moment", "pelvis_list_moment", "pelvis_rotation_moment", "pelvis_tx_force", "pelvis_ty_force", "pelvis_tz_force", "hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment", "knee_angle_moment", "knee_angle_beta_force", "arm_flex_moment", "arm_add_moment", "arm_rot_moment", "ankle_angle_moment", "elbow_flex_moment", "subtalar_angle_moment", "pro_sup_moment", "mtp_angle_moment", "wrist_flex_moment", "wrist_dev_moment"]
        self.results_headers["so"] = ["time", "addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx", "bflh", "bfsh", "edl", "ehl", "fdl", "fhl", "gaslat", "gasmed", "glmax1", "glmax2", "glmax3", "glmed1", "glmed2", "glmed3", "glmin1", "glmin2", "glmin3", "grac", "iliacus", "perbrev", "perlong", "piri", "psoas", "recfem", "sart", "semimem", "semiten", "soleus", "tfl", "tibant", "tibpost", "vasint", "vaslat", "vasmed", "lumbar_ext", "lumbar_bend", "lumbar_rot", "shoulder_flex", "shoulder_add", "shoulder_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["rra"] = []
        self.results_headers["cmc"] = []
        self.results_headers["jr"] = ["time", "hip_on_femur_in_femur_fx", "hip_on_femur_in_femur_fy", "hip_on_femur_in_femur_fz", "hip_on_femur_in_femur_mx", "hip_on_femur_in_femur_my", "hip_on_femur_in_femur_mz", "walker_knee_on_tibia_in_tibia_fx", "walker_knee_on_tibia_in_tibia_fy", "walker_knee_on_tibia_in_tibia_fz", "walker_knee_on_tibia_in_tibia_mx", "walker_knee_on_tibia_in_tibia_my", "walker_knee_on_tibia_in_tibia_mz", "patellofemoral_on_patella_in_patella_fx", "patellofemoral_on_patella_in_patella_fy", "patellofemoral_on_patella_in_patella_fz", "patellofemoral_on_patella_in_patella_mx", "patellofemoral_on_patella_in_patella_my", "patellofemoral_on_patella_in_patella_mz", "ankle_on_talus_in_talus_fx", "ankle_on_talus_in_talus_fy", "ankle_on_talus_in_talus_fz", "ankle_on_talus_in_talus_mx", "ankle_on_talus_in_talus_my", "ankle_on_talus_in_talus_mz", "subtalar_on_calcn_in_calcn_fx", "subtalar_on_calcn_in_calcn_fy", "subtalar_on_calcn_in_calcn_fz", "subtalar_on_calcn_in_calcn_mx", "subtalar_on_calcn_in_calcn_my", "subtalar_on_calcn_in_calcn_mz", "mtp_on_toes_in_toes_fx", "mtp_on_toes_in_toes_fy", "mtp_on_toes_in_toes_fz", "mtp_on_toes_in_toes_mx", "mtp_on_toes_in_toes_my", "mtp_on_toes_in_toes_mz", "back_on_torso_in_torso_fx", "back_on_torso_in_torso_fy", "back_on_torso_in_torso_fz", "back_on_torso_in_torso_mx", "back_on_torso_in_torso_my", "back_on_torso_in_torso_mz"]          
        self.results_headers["bk"] = ["time", "pelvis_X", "pelvis_Y", "pelvis_Z", "pelvis_Ox", "pelvis_Oy", "pelvis_Oz", "femur_X", "femur_Y", "femur_Z", "femur_Ox", "femur_Oy", "femur_Oz", "tibia_X", "tibia_Y", "tibia_Z", "tibia_Ox", "tibia_Oy", "tibia_Oz", "patella_X", "patella_Y", "patella_Z", "patella_Ox", "patella_Oy", "patella_Oz", "talus_X", "talus_Y", "talus_Z", "talus_Ox", "talus_Oy", "talus_Oz", "calcn_X", "calcn_Y", "calcn_Z", "calcn_Ox", "calcn_Oy", "calcn_Oz", "toes_X", "toes_Y", "toes_Z", "toes_Ox", "toes_Oy", "toes_Oz", "torso_X", "torso_Y", "torso_Z", "torso_Ox", "torso_Oy", "torso_Oz", "humerus_X", "humerus_Y", "humerus_Z", "humerus_Ox", "humerus_Oy", "humerus_Oz", "ulna_X", "ulna_Y", "ulna_Z", "ulna_Ox", "ulna_Oy", "ulna_Oz", "radius_X", "radius_Y", "radius_Z", "radius_Ox", "radius_Oy", "radius_Oz", "hand_X", "hand_Y", "hand_Z", "hand_Ox", "hand_Oy", "hand_Oz", "center_of_mass_X", "center_of_mass_Y", "center_of_mass_Z"]


        
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

        # Output database
        self.outfolder = "outputdatabase"
        
        # project
        self.project = "FORCE_SDP"

        # export data
        self.csvfolder = "csvfolder"
        self.csvfileprefix = "force_sdp_results_all_trials"
        self.csvdescfileprefix = "force_sdp_results_subject_descriptives"
        
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

        # additional participant info
        self.additional_participant_info_file = r"outputdatabase\csvfolder\FORCE-ParticipantData-All.csv"
        
        
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
        
        # OpenSim BK parameters
        self.bk_bodies = ["all"]
        self.bk_output_com = True
        self.bk_output_in_local_frame = False


        # ******************************
        # OPENSIM RESULTS
        
        # left leg flip columns (incl. time)
        self.results_flip = {}
        self.results_flip["ik"] = [2, 3, 6, 24, 25]
        self.results_flip["id"] = [2, 3, 6, 14, 15]
        self.results_flip["so"] = []
        self.results_flip["rra"] = []
        self.results_flip["cmc"] = []
        self.results_flip["jr"] = [66, 67, 68, 75, 76, 77, 84, 85, 86, 93, 94, 95, 102, 103, 104, 111, 112, 113, 120, 121, 122]  
        
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
        self.results_columns["jr"] = [[0, 10, 11, 12, 13, 14, 15, 19, 20, 21, 22, 23, 24, 28, 29, 30, 31, 32, 33, 37, 38, 39, 40, 41, 42, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 59, 60, 118, 119, 120, 121, 122, 123],
                                      [0, 64, 65, 66, 69, 68, 69, 73, 74, 75, 76, 77, 78, 82, 83, 84, 85, 86, 87, 91, 92, 93, 94, 95, 96, 100, 101, 102, 103, 104, 105, 109, 110, 111, 112, 113, 114, 118, 119, 120, 121, 122, 123]]
        
        # headers
        self.results_headers = {}
        self.results_headers["ik"] = ["time", "pelvis_tilt", "pelvis_list", "pelvis_rotation", "pelvis_tx", "pelvis_ty", "pelvis_tz", "hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "knee_angle_beta", "ankle_angle", "subtalar_angle", "mtp_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex", "arm_add", "arm_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["id"] = ["time", "pelvis_tilt_moment", "pelvis_list_moment", "pelvis_rotation_moment", "pelvis_tx_force", "pelvis_ty_force", "pelvis_tz_force", "hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment", "knee_angle_moment", "knee_angle_beta_force", "arm_flex_moment", "arm_add_moment", "arm_rot_moment", "ankle_angle_moment", "elbow_flex_moment", "subtalar_angle_moment", "pro_sup_moment", "mtp_angle_moment", "wrist_flex_moment", "wrist_dev_moment"]
        self.results_headers["so"] = ["time", "addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx", "bflh", "bfsh", "edl", "ehl", "fdl", "fhl", "gaslat", "gasmed", "glmax1", "glmax2", "glmax3", "glmed1", "glmed2", "glmed3", "glmin1", "glmin2", "glmin3", "grac", "iliacus", "perbrev", "perlong", "piri", "psoas", "recfem", "sart", "semimem", "semiten", "soleus", "tfl", "tibant", "tibpost", "vasint", "vaslat", "vasmed", "lumbar_ext", "lumbar_bend", "lumbar_rot", "shoulder_flex", "shoulder_add", "shoulder_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        self.results_headers["rra"] = []
        self.results_headers["cmc"] = []
        self.results_headers["jr"] = ["time", "hip_on_femur_in_femur_fx", "hip_on_femur_in_femur_fy", "hip_on_femur_in_femur_fz", "hip_on_femur_in_femur_mx", "hip_on_femur_in_femur_my", "hip_on_femur_in_femur_mz", "walker_knee_on_tibia_in_tibia_fx", "walker_knee_on_tibia_in_tibia_fy", "walker_knee_on_tibia_in_tibia_fz", "walker_knee_on_tibia_in_tibia_mx", "walker_knee_on_tibia_in_tibia_my", "walker_knee_on_tibia_in_tibia_mz", "patellofemoral_on_patella_in_patella_fx", "patellofemoral_on_patella_in_patella_fy", "patellofemoral_on_patella_in_patella_fz", "patellofemoral_on_patella_in_patella_mx", "patellofemoral_on_patella_in_patella_my", "patellofemoral_on_patella_in_patella_mz", "ankle_on_talus_in_talus_fx", "ankle_on_talus_in_talus_fy", "ankle_on_talus_in_talus_fz", "ankle_on_talus_in_talus_mx", "ankle_on_talus_in_talus_my", "ankle_on_talus_in_talus_mz", "subtalar_on_calcn_in_calcn_fx", "subtalar_on_calcn_in_calcn_fy", "subtalar_on_calcn_in_calcn_fz", "subtalar_on_calcn_in_calcn_mx", "subtalar_on_calcn_in_calcn_my", "subtalar_on_calcn_in_calcn_mz", "mtp_on_toes_in_toes_fx", "mtp_on_toes_in_toes_fy", "mtp_on_toes_in_toes_fz", "mtp_on_toes_in_toes_mx", "mtp_on_toes_in_toes_my", "mtp_on_toes_in_toes_mz", "back_on_torso_in_torso_fx", "back_on_torso_in_torso_fy", "back_on_torso_in_torso_fz", "back_on_torso_in_torso_mx", "back_on_torso_in_torso_my", "back_on_torso_in_torso_mz"]          
        
