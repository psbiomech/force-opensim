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
        self.refmodelfile = "LASEM_FORCE_ReferenceModel_Unclamped.osim"

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
        self.subjgroups = ["SYMP", "CTRL"]
        
        # static trial info
        self.staticprefix = "STATIC"
        self.staticused = "Static01"
        self.staticfpchannel = "Fz1"
        
        # file suffixes based on task
        self.trialprefixes = {}
        self.trialprefixes["sdp"] = ["SDP"]
                        
        # file name format regex pattern:
        #   (subjprefix)(num code)_(trialprefix)(num code)
        self.fnpat = "(FAILTCRT|FAILT)(\d+)_([A-Za-z]+)(\d+)([A-Za-z]*)"
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
        