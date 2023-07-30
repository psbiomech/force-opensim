# -*- coding: utf-8 -*-
"""
Run OpenSim pipeline, write input data files

@author: Prasanna Sritharan
"""

import opensim
import scipy.signal as signal
import pandas as pd
import numpy as np
import pickle as pk
import os
import shutil
import re


'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''

# NO CLASSES




'''
-----------------------------------
---- FUNCTIONS: OPENSIM TOOLS -----
-----------------------------------
'''


'''
opensim_pipeline(meta, user):
    Run OpenSim tools pipeline. All data in meta is processed unless specified
    by the restart flag which may have types:
        string: Start from this participant and process until the end
        2-tuple: Process between the first and last participant. To process
                only one participant, set the tuple elements to be the same,
                e.g. ("TRAIL004", "TRAIL004")
'''
def opensim_pipeline(meta, user, analyses, restart = -1):

    # note: I haven't worked out how to run OpenSim in a local folder, simply 
    # changing the pwd doesn't work. So need to clear the log file before each
    # analysis, then copy the log file to the local folder after completion.
    
    # lower case analyses list
    analyses = [a.casefold() for a in analyses]
    
    # run OpenSim for all valid trials
    failedfiles = []
    startflag = 0
    for subj in meta:
        
        # Skip to restart participant, process until last restart participant.
        # Python uses lazy evaluation so combined expressions are efficient.
        if restart != -1:
            if startflag == 1:
                if (type(restart) == tuple) and (subj == restart[1]):
                    startflag = 0            
            elif startflag == 0:
                if (type(restart) == str) and (subj == restart):
                    startflag = 1
                elif (type(restart) == tuple) and (subj == restart[0]):
                    if restart[0] != restart[1]:
                        startflag = 1
                else:
                    continue        
        
        for group in meta[subj]["trials"]:
            
            # copy analyses list for subject
            analyses0 = [ans.casefold() for ans in analyses]
            
            # find the static trial for model scaling
            for trial in meta[subj]["trials"][group]:
                            
                #****** FOR TESTING ONLY ******
                # trialre = re.compile("SKIPPING THIS!")
                # if not trialre.match(trial):
                #     # print("%s ---> SKIP" % trial)
                #     continue
                #******************************

                # run the scale tool if the static trial is the trial to be
                # used for scaling
                if "scale" in analyses0:
                    if meta[subj]["trials"][group][trial]["usedstatic"]:
    
                        try:
                        
                            # load the OsimKey
                            pklpath = meta[subj]["trials"][group][trial]["outpath"]
                            pklfile = meta[subj]["trials"][group][trial]["trial"] + "_osimkey.pkl"
                            with open(os.path.join(pklpath, pklfile),"rb") as fid: 
                                osimkey = pk.load(fid)
                                
                            # create an OpenSim log folder
                            logfolder = os.path.join(pklpath, user.triallogfolder)
                            if not os.path.isdir(logfolder):
                                os.makedirs(logfolder)
                            
                            # run the scale tool
                            run_opensim_scale(osimkey, user)
                            
                            # get the full model path
                            modelfile = meta[subj]["trials"][group][trial]["osim"]
                            modelfullpathfile = os.path.join(pklpath, modelfile)
                            
                            # copy the model into the subject root folder
                            subjrootpath = os.path.join(user.rootpath, user.outfolder, subj)
                            shutil.copy(modelfullpathfile, subjrootpath)                            
                            
                            # once a model is created we can stop finding valid
                            # static trials
                            analyses0.remove("scale")
                            break
                
                        except:
                            print("%s ---> ***FAILED***" % trial)
                            failedfiles.append(trial)
                            #raise
            
            # if scale was the only analysis, then go to next group
            if not analyses0: continue
            
            # find dynamic trials, always copy model into folder and run the
            # requested analyses
            for trial in meta[subj]["trials"][group]:
                
                #****** FOR TESTING ONLY ******
                # trialre = re.compile("FAILTCRT01_SLDJ01")
                # if not trialre.match(trial):
                #     #print("%s ---> SKIP" % trial)
                #     continue
                #******************************
                
                if not meta[subj]["trials"][group][trial]["isstatic"]:

                    try:
                        
                        # load the OsimKey
                        pklpath = meta[subj]["trials"][group][trial]["outpath"]
                        pklfile = meta[subj]["trials"][group][trial]["trial"] + "_osimkey.pkl"
                        with open(os.path.join(pklpath, pklfile),"rb") as fid: 
                            osimkey = pk.load(fid)
    
                        # create an OpenSim log folder
                        logfolder = os.path.join(pklpath, user.triallogfolder)
                        if not os.path.isdir(logfolder):
                            os.makedirs(logfolder)
                        
                        # copy the model into the trial folder
                        modelfile = meta[subj]["trials"][group][trial]["osim"]
                        subjrootpathfile = os.path.join(user.rootpath, user.outfolder, subj, modelfile)
                        shutil.copy(subjrootpathfile, pklpath)
                        
                        # run the required analyses
                        for ans in analyses0:
                            
                            # skip scale
                            if ans == "scale": continue
                            
                            # create output folder
                            if not os.path.exists(os.path.join(pklpath, ans)): os.makedirs(os.path.join(pklpath, ans))
                            
                            # analyses
                            if ans == "ik":
                                run_opensim_ik(osimkey, user)
                            elif ans == "id":
                                run_opensim_id(osimkey, user)
                            elif ans == "so":
                                run_opensim_so(osimkey, user)
                            elif ans == "rra":
                                run_opensim_rra(osimkey, user)
                            elif ans == "cmc":
                                run_opensim_cmc(osimkey, user)
                            elif ans == "jr":
                                run_opensim_jr(osimkey, user)
                            elif ans == "bk":
                                    run_opensim_bk(osimkey, user)
                                
                    except:
                        print("%s ---> ***FAILED***" % trial)
                        failedfiles.append(trial)
                        #raise
                            
    return failedfiles



'''
run_opensim_scale(osimkey, user):    
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run 
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_scale(osimkey, user):
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    subject = osimkey.subject
    model = osimkey.model
    trial = osimkey.trial

    # clear log file
    open(user.logfile, "w").close()

    print("\nCreating scaled model: %s" % subject)
    print("------------------------------------------------")
    
    # create an ScaleTool from a generic setup file
    print("Create new ScaleTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupscale
    tool = opensim.ScaleTool(os.path.join(refsetuppath, refsetupfile))
    tool.setPathToSubject("")
    
    # set subject mass
    tool.setSubjectMass(osimkey.mass)
    
    # set subject name
    tool.setName(osimkey.subject)
    
    
    # ******************************
    # GENERIC MODEL MAKER
    
    print("Initialising GenericModelMaker...")
    
    # set the model file name
    refmodelpath = user.refmodelpath
    refmodelfile = user.refmodelfile
    modelmaker = tool.getGenericModelMaker()
    modelmaker.setModelFileName(os.path.join(refmodelpath, refmodelfile))
    

    # ******************************
    # MODEL SCALER
    
    print("Setting up ModelScaler...")
    
    # set the static TRC file to be used for scaling
    modelscaler = tool.getModelScaler()
    modelscaler.setMarkerFileName(os.path.join(fpath, trial + "_markers.trc"))
    
    # set time window for scaling to 45%-55% of trial time
    t0 = osimkey.events["time"][0] + (osimkey.events["time"][1] - osimkey.events["time"][0]) * 0.45
    t1 = osimkey.events["time"][0] + (osimkey.events["time"][1] - osimkey.events["time"][0]) * 0.55
    twindow = opensim.ArrayDouble(0, 2)
    twindow.set(0, t0)
    twindow.set(1, t1)
    modelscaler.setTimeRange(twindow)
    
    # set output model file name
    modelscaler.setOutputModelFileName(os.path.join(fpath, model))
    
    
    # ******************************
    # MARKER PLACER   

    print("Setting up MarkerPlacer...")

    # set the static TRC file to be used for scaling
    markerplacer = tool.getMarkerPlacer()
    markerplacer.setMarkerFileName(os.path.join(fpath, trial + "_markers.trc"))
    
    # set time window
    markerplacer.setTimeRange(twindow)
    
    # set static trial output motion file
    markerplacer.setOutputMotionFileName(os.path.join(fpath, trial + "_static_ik.mot"))
    
    # set output model file
    markerplacer.setOutputModelFileName(os.path.join(fpath, model))    
     
    # set output marker file
    markerplacer.setOutputMarkerFileName(os.path.join(fpath, trial + "_markers.xml"))   


    # ******************************
    # RUN TOOL
    
    print("Running the ScaleTool...")
    
    # save the settings in a setup file
    tool.printToXML(os.path.join(fpath, trial + "_Setup_Scale.xml"))

    # run the tool, and scale the model if required
    try:
        
        # run the tool
        tool.run()
        
        
        # ******************************
        # UPDATE MODEL MUSCLE PROPERTIES
        # (comment out any unused code blocks below as required)
        
        # scale the model FoM if required
        if hasattr(user, "fom_scalefactor"):
            sf_fom = user.fom_scalefactor
            if (type(sf_fom) is dict):
                print("---> Scaling muscle FoM in model...")
                shutil.copyfile(os.path.join(fpath, model), os.path.join(fpath, subject + "_original_FoM.osim"))
                model0 = opensim.Model(os.path.join(fpath, model))
                refmodel = opensim.Model(os.path.join(refmodelpath, refmodelfile))
                model1 = update_osim_fom(model0, sf_fom, refmodel)
                model1.printToXML(os.path.join(fpath, model))
            
        # scale the model LoM (maintain constant LMT) if required
        if hasattr(user, "lom_lmt_scalefactor"):
            sf_lom_lmt = user.lom_lmt_scalefactor
            if (type(sf_lom_lmt) is dict):
                print("---> Scaling muscle LoM in model (LMT remains constant)...")
                shutil.copyfile(os.path.join(fpath, model), os.path.join(fpath, subject + "_original_LoM_LsT.osim"))
                model0 = opensim.Model(os.path.join(fpath, model))
                model1 = update_osim_lom_const_lmt(model0, sf_lom_lmt)
                model1.printToXML(os.path.join(fpath, model))  

        # scale the model LoM if required
        if hasattr(user, "lom_scalefactor"):
            sf_lom = user.lom_scalefactor
            if (type(sf_lom) is dict):
                print("---> Scaling muscle LoM in model...")
                shutil.copyfile(os.path.join(fpath, model), os.path.join(fpath, subject + "_original_LoM.osim"))
                model0 = opensim.Model(os.path.join(fpath, model))
                model1 = update_osim_lom(model0, sf_lom)
                model1.printToXML(os.path.join(fpath, model)) 

        # scale the model LsT if required
        if hasattr(user, "lst_scalefactor"):
            sf_lst = user.lst_scalefactor
            if (type(sf_lst) is dict):
                print("---> Scaling muscle LsT in model...")
                shutil.copyfile(os.path.join(fpath, model), os.path.join(fpath, subject + "_original_LsT.osim"))
                model0 = opensim.Model(os.path.join(fpath, model))
                model1 = update_osim_lom_const_lmt(model0, sf_lst)
                model1.printToXML(os.path.join(fpath, model))   


        print("Done.")
                
    except:
        print("---> ERROR: Scale failed. Skipping Scale for %s." % trial)
    finally:
        shutil.copyfile(user.logfile, os.path.join(fpath, user.triallogfolder, "out_SCALE.log"))
        print("------------------------------------------------\n")
    
    # ******************************    
    
    return None



'''
run_opensim_ik(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_ik(osimkey, user):

    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    # clear log file
    open(user.logfile, "w").close()

    print("\nPerforming IK on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an IK Tool from a generic setup file
    print("Create new IKTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupik
    tool = opensim.InverseKinematicsTool(os.path.join(refsetuppath, refsetupfile))
    
    # load the model
    print("Loading the model: %s..." % modelfile)
    model = opensim.Model(os.path.join(fpath, modelfile))
    model.initSystem()
    
    # set the model in the tool
    tool.setModel(model)

    # set the initial and final times (limit to between first and last event,
    # but allow extra time at start and end for RRA/CMC)
    t0 = float(osimkey.events["time"][0]) + user.rra_start_time_offset
    t1 = float(osimkey.events["time"][-1]) + user.rra_end_time_offset
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setStartTime(t0)
    tool.setEndTime(t1)

    # set input TRC file
    print("Setting the marker data (TRC) file...")
    tool.setMarkerDataFileName(os.path.join(fpath, trial + "_markers.trc"))
    
    # output MOT file location
    print("Setting output file name...")
    motfilepath = os.path.join(fpath, user.ikcode)
    if not os.path.isdir(motfilepath): os.makedirs(motfilepath)
    tool.setOutputMotionFileName(os.path.join(motfilepath, trial + "_ik.mot"))  


    # ******************************
    # RUN TOOL

    print("Running the IKTool...")

    # save the settings in a setup file
    tool.printToXML(os.path.join(fpath, trial + "_Setup_IK.xml"))
        
    # run the tool
    try:
        tool.run()     
        print("Done.")
    except:
        print("---> ERROR: IK failed. Skipping IK for %s." % trial)
    finally:
        shutil.copyfile(user.logfile, os.path.join(fpath, user.triallogfolder, "out_IK.log")) 
        print("------------------------------------------------\n")
    
    # ******************************
    
    return None



'''
run_opensim_id(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_id(osimkey, user):
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    # clear log file
    open(user.logfile, "w").close()
    
    print("\nPerforming ID on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an ID Tool from a generic setup file
    print("Create new IDTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupid
    tool = opensim.InverseDynamicsTool(os.path.join(refsetuppath, refsetupfile))

    # load the model
    print("Loading the model: %s..." % modelfile)
    model = opensim.Model(os.path.join(fpath, modelfile))
    model.initSystem()
    
    # set the model in the tool
    tool.setModel(model)   
    
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0])
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]])
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setStartTime(t0)
    tool.setEndTime(t1)

    # set input directory and coordinates data file
    print("Setting coordinates data file...")
    tool.setCoordinatesFileName(os.path.join(fpath, user.ikcode, trial + "_ik.mot"))
    tool.setLowpassCutoffFrequency(6.0)
    
    # set output directories and generalised forces storage file (note:
    # InverseDynamicsTool XML parser does not seem to like full paths
    # for the OutputGenForceFileName tag, so need to set results dir)
    print("Setting output file name...")
    stofilepath = os.path.join(fpath, user.idcode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)
    tool.setOutputGenForceFileName(trial + "_id.sto")
    
    # create an external loads object from template, set it up, and
    # print to destination folder (This is not the best way to do this,
    # but Opensim Tools are designed to read directly from XML files,
    # so better to fully set up an external loads file, print it then
    # load it again into the Tool, than to create an ExternalLoads
    # object and connect it to the Model. This also ensures a copy of
    # the external loads file is available in the trial folder in case
    # a one-off analysis needs to be run in future.)
    print("Creating external loads XML file...")
    extloadsfile = os.path.join(fpath, trial + "_ExternalLoads.xml")
    if not os.path.isfile(extloadsfile):
        extloads = opensim.ExternalLoads(os.path.join(refsetuppath, user.additionalfilesfolder, user.refexternalloads), True)       
        extloads.setDataFileName(os.path.join(fpath, trial + "_grf.mot"))
        extloads.printToXML(extloadsfile)  

    # set the external loads file name
    tool.setExternalLoadsFileName(extloadsfile)
    
    
    # ******************************
    # RUN TOOL 
    
    print("Running the IDTool...")

    # save the settings in a setup file
    tool.printToXML(os.path.join(fpath, trial + "_Setup_ID.xml"))
    
    # run the tool
    try:
        tool.run()
        print("Done.")
    except:
        print("---> ERROR: ID failed. Skipping ID for %s." % trial)
    finally:
        shutil.copyfile(user.logfile, os.path.join(fpath, user.triallogfolder, "out_ID.log")) 
        print("------------------------------------------------\n")

    # ******************************
        
    return None    
   

    
'''
run_opensim_bk(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_bk(osimkey, user):
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    # clear log file
    open(user.logfile, "w").close()
    
    print("\nPerforming BK on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an generic AnalyzeTool from the template setup file (set the flag
    # aLoadModelAndInput = false as we are only configuring the setup file at 
    # this stage)
    print("Create new AnalyzeTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupso
    tool = opensim.AnalyzeTool(os.path.join(refsetuppath, refsetupfile), False)
    tool.setName(trial) 
  
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0])
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]])
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)
    
    # set output directory
    print("Setting output file name...")
    stofilepath = os.path.join(fpath, user.bkcode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)
    
    # # create an external loads object from template, set it up, and
    # # print to destination folder (This is not the best way to do this,
    # # but Opensim Tools are designed to read directly from XML files,
    # # so better to fully set up an external loads file, print it then
    # # load it again into the Tool, than to create an ExternalLoads
    # # object and connect it to the Model. This also ensures a copy of
    # # the external loads file is available in the trial folder in case
    # # a one-off analysis needs to be run in future.)
    # print("Creating external loads XML file...")
    # extloadsfile = os.path.join(fpath, trial + "_ExternalLoads.xml")
    # if not os.path.isfile(extloadsfile):
    #     extloads = opensim.ExternalLoads(os.path.join(refsetuppath, user.additionalfilesfolder, user.refexternalloads), True)       
    #     extloads.setDataFileName(os.path.join(fpath, trial + "_grf.mot"))
    #     extloads.printToXML(extloadsfile)  

    # # set the external loads file name
    # tool.setExternalLoadsFileName(extloadsfile)
    
    
    # ******************************
    # CREATE BODY KINEMATICS ANALYSIS   
    
    # create an JR analysis
    print("Create new BK Analysis to the AnalysisSet...")
    bk = opensim.BodyKinematics()
    bk.setName("bk")
    bk.setStartTime(t0)
    bk.setEndTime(t1)
    
    # # set the muscle forces 
    # print("Setting muscle forces...")
    # if user.jr_use_cmc_forces:
    #     forcefile = os.path.join(fpath, user.cmccode, trial + "_Actuation_force.sto")
    # else:
    #     forcefile = os.path.join(fpath, user.socode, trial + "_so_force.sto")         
    # jr.setForcesFileName(forcefile)
    
    # set the bodies to be analysed
    print("Setting bodies to be analysed...")
    bodies = opensim.ArrayStr()
    for j in user.bk_bodies:
        bodies.append(j)
    bk.setBodiesToRecord(bodies)
    bk.setRecordCenterOfMass(user.bk_output_com)

    # Set frame of reference
    bk.setExpressResultsInLocalFrame(user.bk_output_in_local_frame)    


    # add the BodyKinematics analysis to the tool AnalysisSet
    print("Append new BK Analysis to the AnalysisSet...")
    analyses = tool.getAnalysisSet()
    analyses.insert(0, bk)


    # ******************************
    # SET THE MODEL AND KINEMATICS
    # (use baseline or RRA adjusted model and kinematics)
           
    # desired model name
    print("Loading the model...")
    if user.bk_use_cmc_results:
        if user.update_mass:
            actualmodelfile = trial + "_RRA_" + str(user.rraiter) + "_TorsoAdjusted_MassUpdated.osim"
        else:
            actualmodelfile = trial + "_RRA_" + str(user.rraiter) + "_TorsoAdjusted.osim"
    else:
        actualmodelfile = modelfile
    tool.setModelFilename(os.path.join(fpath, actualmodelfile))
        
    # set coordinates data file
    print("Setting coordinates data file...")
    if user.jr_use_cmc_forces:
        kinfile = os.path.join(fpath, user.cmccode, trial + "_Kinematics_q.sto")
        filtfreq = -1    
    else:            
        kinfile = os.path.join(fpath, user.ikcode, trial + "_ik.mot")
        filtfreq = 6.0
    tool.setCoordinatesFileName(kinfile)
    tool.setLowpassCutoffFrequency(filtfreq)
    
    
    # ******************************
    # RUN TOOL 
    
    print("Running the AnalysisTool (BK)...")

    # save the settings in a setup file
    customsetupfile = os.path.join(fpath, trial + "_Setup_BK.xml")
    tool.printToXML(customsetupfile)
    
    # run the tool (need to load the setup again into a new AnalyzeTool)
    try:
        tool2 = opensim.AnalyzeTool(customsetupfile)
        tool2.run()
        print("Done.")
    except:
        print("---> ERROR: BK failed. Skipping BK for %s." % trial)
    finally:
        shutil.copyfile(user.logfile, os.path.join(fpath, user.triallogfolder, "out_BK.log")) 
        print("------------------------------------------------\n")

    # ******************************
        
    return None   



'''
run_opensim_so(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_so(osimkey, user):
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    # clear log file
    open(user.logfile, "w").close()
    
    print("\nPerforming SO on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an generic AnalyzeTool from the template setup file (set the flag
    # aLoadModelAndInput = false as we are only configuring the setup file at 
    # this stage)
    print("Create new AnalyzeTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupso
    tool = opensim.AnalyzeTool(os.path.join(refsetuppath, refsetupfile), False)
    tool.setName(trial)
    
    # set the model in the tool
    print("Loading the model: %s..." % modelfile)
    tool.setModelFilename(os.path.join(fpath, modelfile))   
  
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0])
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]])
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)

    # set coordinates data file
    print("Setting coordinates data file...")
    tool.setCoordinatesFileName(os.path.join(fpath, user.ikcode, trial + "_ik.mot"))
    tool.setLowpassCutoffFrequency(6.0)
    
    # set output directory
    print("Setting output file name...")
    stofilepath = os.path.join(fpath, user.socode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)
    
    # create an external loads object from template, set it up, and
    # print to destination folder (This is not the best way to do this,
    # but Opensim Tools are designed to read directly from XML files,
    # so better to fully set up an external loads file, print it then
    # load it again into the Tool, than to create an ExternalLoads
    # object and connect it to the Model. This also ensures a copy of
    # the external loads file is available in the trial folder in case
    # a one-off analysis needs to be run in future.)
    print("Creating external loads XML file...")
    extloadsfile = os.path.join(fpath, trial + "_ExternalLoads.xml")
    if not os.path.isfile(extloadsfile):
        extloads = opensim.ExternalLoads(os.path.join(refsetuppath, user.additionalfilesfolder, user.refexternalloads), True)       
        extloads.setDataFileName(os.path.join(fpath, trial + "_grf.mot"))
        extloads.printToXML(extloadsfile)  

    # set the external loads file name
    tool.setExternalLoadsFileName(extloadsfile)
    
    
    # ******************************
    # CREATE STATIC OPTIMISATION ANALYSIS   
    
    # create an SO analysis
    print("Append new SO Analysis to the AnalysisSet...")
    so = opensim.StaticOptimization()
    so.setName("so")
    so.setStartTime(t0)
    so.setEndTime(t1)
    
    # add the StaticOptimization analysis to the tool AnalysisSet
    analyses = tool.getAnalysisSet()
    analyses.insert(0, so)


    # ******************************
    # PREPARE RESERVE ACTUATORS
    # (set the pelvis reserves to act at pelvis COM, or use default application
    # at the ground-pelvis origin located at pelvis geometric centre)
        
    # get the pelvis COM location from the model
    model = opensim.Model(os.path.join(fpath, modelfile))
    pelvis = opensim.Body.safeDownCast(model.findComponent("pelvis"))
    pelviscom = pelvis.getMassCenter() 
    
    # load reference actuator forceset, get the pelvis actuators and set force
    # application point to pelvis COM
    refrefreserveactuators = user.refreserveactuators
    residforceset = opensim.ForceSet(os.path.join(refsetuppath, user.additionalfilesfolder, refrefreserveactuators))
    for x in ["FX","FY","FZ"]:
        residforce = opensim.PointActuator.safeDownCast(residforceset.get(x))
        residforce.set_point(pelviscom)

    # write updated actuator set to file    
    forcesetfile = trial + "_Reserve_Actuators.xml"
    residforceset.printToXML(os.path.join(fpath, forcesetfile))
      
    # append reserve actuators
    # (note: need to create a new ArrayStr with one element, add the file name,
    # and then pass the ArrayStr to setForceSetFiles; however, the tool does
    # not like whitespace in the <force_set_files> XML property)
    print("Appending reserve actuators...")
    fsvec = opensim.ArrayStr()    
    fsvec.append(forcesetfile)
    tool.setForceSetFiles(fsvec)
    tool.setReplaceForceSet(False)   

    # # append default reference actuators, with pelvis residuals acting at 
    # # ground-pelvis origin located at pelvis geometric centre
    # refforcesetfile = user.refreserveactuators
    # residforceset = opensim.ForceSet(os.path.join(refsetuppath, refforcesetfile))
    # forcesetfile = trial + "_Reserve_Actuators.xml"
    # residforceset.printToXML(os.path.join(fpath, forcesetfile))
    # fsvec = opensim.ArrayStr()
    # fsvec.append(forcesetfile)
    # tool.setForceSetFiles(fsvec)
    # tool.setReplaceForceSet(False)     

    
    # ******************************
    # RUN TOOL 
    
    print("Running the AnalysisTool (SO)...")

    # save the settings in a setup file
    customsetupfile = os.path.join(fpath, trial + "_Setup_SO.xml")
    tool.printToXML(customsetupfile)
    
    # run the tool (need to load the setup again into a new AnalyzeTool)
    try:
        tool2 = opensim.AnalyzeTool(customsetupfile)
        tool2.run()
        print("Done.")
    except:
        print("---> ERROR: SO failed. Skipping SO for %s." % trial)
    finally:
        shutil.copyfile(user.logfile, os.path.join(fpath, user.triallogfolder, "out_SO.log")) 
        print("------------------------------------------------\n")

    # ******************************
        
    return None    
    


'''
run_opensim_rra(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_rra(osimkey, user):

    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    # clear log file
    open(user.logfile, "w").close()
    
    print("\nPerforming RRA on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an RRA Tool from a generic setup file
    # (note: set loadModel=false for manual load later)
    print("Create new RRATool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetuprra
    tool = opensim.RRATool(os.path.join(refsetuppath, refsetupfile), False)      
    
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0]) + user.rra_start_time_offset
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]]) + user.rra_end_time_offset
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)

    # set flag to adjust model COM and set body to torso
    tool.setAdjustCOMToReduceResiduals(True)
    tool.setAdjustedCOMBody("torso")

    # set desired kinematics file name and filter frequency
    tool.setDesiredKinematicsFileName(os.path.join(fpath, user.ikcode, trial + "_ik.mot"))
    tool.setLowpassCutoffFrequency(user.kinematics_filter_cutoff)
    
    # create an external loads object from template, set it up, and
    # print to destination folder (This is not the best way to do this,
    # but Opensim Tools are designed to read directly from XML files,
    # so better to fully set up an external loads file, print it then
    # load it again into the Tool, than to create an ExternalLoads
    # object and connect it to the Model. This also ensures a copy of
    # the external loads file is available in the trial folder in case
    # a one-off analysis needs to be run in future.)
    print("Creating external loads XML file...")
    extloadsfile = os.path.join(fpath, trial + "_ExternalLoads.xml")
    if not os.path.isfile(extloadsfile):
        extloads = opensim.ExternalLoads(os.path.join(refsetuppath, user.additionalfilesfolder, user.refexternalloads), True)       
        extloads.setDataFileName(os.path.join(fpath, trial + "_grf.mot"))
        extloads.printToXML(extloadsfile)  
    
    # set the external loads file name
    tool.setExternalLoadsFileName(extloadsfile)
    
    # set output directory
    print("Setting output folder name...")
    stofilepath = os.path.join(fpath, user.rracode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)


    # ******************************
    # PREPARE RRA ACTUATORS
    # (set the pelvis reserves to act at pelvis COM, or use default application
    # at the ground-pelvis origin located at pelvis geometric centre)
    
    # get the pelvis COM location from the model
    model = opensim.Model(os.path.join(fpath, modelfile))
    pelvis = opensim.Body.safeDownCast(model.findComponent("pelvis"))
    pelviscom = pelvis.getMassCenter()
    
    # load reference actuator forceset, get the pelvis actuators and set force
    # application point to pelvis COM
    refforcesetfile = user.refrraactuators
    residforceset = opensim.ForceSet(os.path.join(refsetuppath, user.additionalfilesfolder, refforcesetfile))
    for x in ["FX","FY","FZ"]:
        residforce = opensim.PointActuator.safeDownCast(residforceset.get(x))
        residforce.set_point(pelviscom)

    # write actuator set to file
    forcesetfile = trial + "_RRA_Actuators.xml"
    residforceset.printToXML(os.path.join(fpath, forcesetfile))
    
    # replace reserve actuators
    # (note: need to create a new ArrayStr with one element, add the file name,
    # and then pass the ArrayStr to setForceSetFiles; however, the tool does
    # not like whitespace in the <force_set_files> XML property)
    fsvec = opensim.ArrayStr();
    fsvec.append(forcesetfile)
    tool.setForceSetFiles(fsvec)
    tool.setReplaceForceSet(True)
    
    # # replace actuators with default reference actuators, with pelvis residuals
    # # acting at ground-pelvis origin located at pelvis geometric centre
    # refforcesetfile = user.refrraactuators
    # residforceset = opensim.ForceSet(os.path.join(refsetuppath, refforcesetfile))
    # forcesetfile = trial + "_RRA_Actuators.xml"
    # residforceset.printToXML(os.path.join(fpath, forcesetfile))    
    # fsvec = opensim.ArrayStr()
    # fsvec.append(forcesetfile)
    # tool.setForceSetFiles(fsvec)
    # tool.setReplaceForceSet(True)    
    
    
    # ******************************
    # PREPARE RRA TASKS
 
    # load RRA tasks set
    rratasksfile = user.refrratasks
    rrataskset = opensim.CMC_TaskSet(os.path.join(refsetuppath, user.additionalfilesfolder, rratasksfile))
    
    # updates here...
    # e.g. controller weights and gains
    
    # print tasks to trial folder
    rrataskset.printToXML(os.path.join(fpath, trial + "_RRA_Tasks.xml"))
    
    # set RRA tasks file in tool
    tool.setTaskSetFileName(os.path.join(fpath, trial + "_RRA_Tasks.xml"))
       
    
    # ******************************
    # RUN TOOL
    
    print("Running the RRATool, %d iterations...\n---> (this may take a while)" %  user.rraiter)

    try:
        
        # run the tool require number of iterations
        rra_adjusted_model_file = []
        for i in range(user.rraiter):
            
            # clear log file
            open(user.logfile, "w").close()
            
            # set tool name based on current iteration
            rraiter = "RRA_" + str(i + 1)
            rraname = trial + "_" + rraiter
            tool.setName(rraname)
            
            # get the model, ensure model MTP joints are locked (only needs to
            # be done on the first iteration)
            if i==0:
                rramodelfile = os.path.join(fpath, modelfile)
            else:
                rramodelfile = rra_adjusted_model_file 

            # update the model to prescribe upper body kinematics from IK
            if (i == 0) and user.prescribe_upper_body_motion:
                coord_list = user.prescribe_coord_list
                model0 = opensim.Model(rramodelfile)   
                ikdata = pd.read_csv(os.path.join(fpath, user.ikcode, trial + "_ik.mot"), sep = "\t", header = 8)
                model1 = prescribe_kinematics(model0, ikdata, coord_list, user.kinematics_filter_butter_order, user.kinematics_filter_cutoff)
                rramodelfile = rramodelfile.rstrip(".osim") + "_PrescribedUpper.osim"
                model1.printToXML(rramodelfile)

            # set model file name
            tool.setModelFilename(rramodelfile)
            
            # set new adjusted model file name
            rra_adjusted_model_file = os.path.join(fpath, rraname + "_TorsoAdjusted.osim")
            tool.setOutputModelFileName(rra_adjusted_model_file)
            
            # save the current settings in a setup file
            rrasetupfile = os.path.join(fpath, trial + "_Setup_" + rraiter + ".xml")
            tool.printToXML(rrasetupfile)
            
            print("---> Iteration: %d" % (i + 1))
            
            # load the current setup file, run the current RRA tool
            try:               
                rratool2 = opensim.RRATool(rrasetupfile)
                rratool2.run()
            except:
                print("---> ERROR: RRA iteration %d failed. Skipping RRA for %s." % (i + 1, trial))
            else:
                shutil.copyfile(user.logfile, os.path.join(fpath, user.triallogfolder, "out_" + rraiter + ".log"))                     

 
            # update the RRA model segment masses, print the new model and
            # overwrite the model file parameter for the next iteration
            if user.update_mass:
                rra_model = opensim.Model(rra_adjusted_model_file)
                rra_adjusted_model = perform_recommended_mass_change(rra_model, os.path.join(fpath, user.triallogfolder, "out_" + rraiter + ".log"))
                rra_adjusted_model_file = os.path.join(fpath, rraname + "_TorsoAdjusted_MassUpdated.osim")
                rra_adjusted_model.printToXML(rra_adjusted_model_file)
            

    except: 
        print("---> ERROR: RRA failed. Skipping RRA for %s. Log file not copied." % trial)
    else:
        print("Done.")
    finally:
        print("------------------------------------------------\n")
    
    # ******************************
        
    return None    
        


'''
run_opensim_cmc(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_cmc(osimkey, user):

    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    # clear log file
    open(user.logfile, "w").close()
    
    print("\nPerforming CMC on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create a CMC Tool from a generic setup file
    # (note: set loadModel=false for manual load later)
    print("Create new CMCTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupcmc
    tool = opensim.CMCTool(os.path.join(refsetuppath, refsetupfile), False)
    tool.setName(trial)     
    
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0]) + user.cmc_start_time_offset
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]]) + user.cmc_end_time_offset
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)
    
    # create an external loads object from template, set it up, and
    # print to destination folder (This is not the best way to do this,
    # but Opensim Tools are designed to read directly from XML files,
    # so better to fully set up an external loads file, print it then
    # load it again into the Tool, than to create an ExternalLoads
    # object and connect it to the Model. This also ensures a copy of
    # the external loads file is available in the trial folder in case
    # a one-off analysis needs to be run in future.)
    print("Creating external loads XML file...")
    extloadsfile = os.path.join(fpath, trial + "_ExternalLoads.xml")
    if not os.path.isfile(extloadsfile):
        extloads = opensim.ExternalLoads(os.path.join(refsetuppath, user.additionalfilesfolder, user.refexternalloads), True)       
        extloads.setDataFileName(os.path.join(fpath, trial + "_grf.mot"))
        extloads.printToXML(extloadsfile)  
    
    # set the external loads file name
    tool.setExternalLoadsFileName(extloadsfile)
    
    # set output directory
    print("Setting output folder name...")
    stofilepath = os.path.join(fpath, user.cmccode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)

    # use fast target (should normally be set to true, but set to false if 
    # having difficulty solving, e.g. pathological movement patterns or dynamic
    # movement patterns such as fast running)
    tool.setUseFastTarget(user.use_fast_target)


    # ******************************
    # PREPARE CMC ACTUATORS
    # (set the pelvis reserves to act at pelvis COM, or use default application
    # at the ground-pelvis origin located at pelvis geometric centre)
    
    # get the pelvis COM location from the desired model
    # (the pelvis COM should be the same for all models anyway)
    if user.use_rra_model:
        rramodelfile = trial + "_RRA_" + str(user.rraiter) + "_TorsoAdjusted.osim"
    else:
        rramodelfile = modelfile   
    model = opensim.Model(os.path.join(fpath, rramodelfile))
    pelvis = opensim.Body.safeDownCast(model.findComponent("pelvis"))
    pelviscom = pelvis.getMassCenter()
    
    # load reference actuator forceset, get the pelvis actuators and set force
    # application point to pelvis COM
    refforcesetfile = user.refcmcactuators
    residforceset = opensim.ForceSet(os.path.join(refsetuppath, user.additionalfilesfolder, refforcesetfile))
    for x in ["FX","FY","FZ"]:
        residforce = opensim.PointActuator.safeDownCast(residforceset.get(x))
        residforce.set_point(pelviscom)

    # write actuator set to file
    print("Preparing CMC actuators...")
    forcesetfile = trial + "_CMC_Actuators.xml"
    residforceset.printToXML(os.path.join(fpath, forcesetfile))
    
    # replace reserve actuators
    # (note: need to create a new ArrayStr with one element, add the file name,
    # and then pass the ArrayStr to setForceSetFiles; however, the tool does
    # not like whitespace in the <force_set_files> XML property)
    fsvec = opensim.ArrayStr();
    fsvec.append(forcesetfile)
    tool.setForceSetFiles(fsvec)
    tool.setReplaceForceSet(False)
    
    # # append default reference actuators, with pelvis residuals acting at 
    # # ground-pelvis origin located at pelvis geometric centre
    # refforcesetfile = user.refcmcactuators
    # residforceset = opensim.ForceSet(os.path.join(refsetuppath, refforcesetfile))
    # forcesetfile = trial + "_CMC_Actuators.xml"
    # residforceset.printToXML(os.path.join(fpath, forcesetfile))      
    # fsvec = opensim.ArrayStr()
    # fsvec.append(forcesetfile)
    # tool.setForceSetFiles(fsvec)
    # tool.setReplaceForceSet(False)    
    
    
    # ******************************
    # PREPARE CMC TASKS
    # (allows for trial-specific controller weights and gains to be applied)
 
    # load CMC tasks set
    print("Preparing CMC tasks...")
    cmctasksfile = user.refcmctasks
    cmctaskset = opensim.CMC_TaskSet(os.path.join(refsetuppath, user.additionalfilesfolder, cmctasksfile))
    
    # updates here... (TBD)
    # e.g. controller weights and gains
    
    # print tasks to trial folder
    cmctaskset.printToXML(os.path.join(fpath, trial + "_CMC_Tasks.xml"))
    
    # set CMC tasks file in tool
    tool.setTaskSetFileName(os.path.join(fpath, trial + "_CMC_Tasks.xml"))
       

    # ******************************
    # PREPARE CMC CONTROL CONSTRAINTS
    # (allows for trial-specific controller constraints to be applied)

    # load reference CMC control constraints
    print("Preparing CMC control constraints...")
    cmccontrolsfile = user.refcmccontrolconstraints
    cmccontrolset = opensim.ControlSet(os.path.join(refsetuppath, user.additionalfilesfolder, cmccontrolsfile))
    
    # updates here... (TBD)
    # e.g. actuator gains
    
    # print control constraints to trial folder
    cmccontrolset.printToXML(os.path.join(fpath, trial + "_CMC_ControlConstraints.xml"))
    
    # set CMC control constraints file in tool
    tool.setConstraintsFileName(os.path.join(fpath, trial + "_CMC_ControlConstraints.xml"))
    

    # ******************************
    # SET THE MODEL AND KINEMATICS
    # (use baseline or RRA adjusted model and kinematics)

    # desired model name
    print("Loading the model...")
    if user.use_rra_model:
        if user.update_mass:
            actualmodelfile = trial + "_RRA_" + str(user.rraiter) + "_TorsoAdjusted_MassUpdated.osim"
        else:
            actualmodelfile = trial + "_RRA_" + str(user.rraiter) + "_TorsoAdjusted.osim"
    else:
        actualmodelfile = modelfile
    tool.setModelFilename(os.path.join(fpath, actualmodelfile))
        
    # set desired kinematics file and filter frequency
    print("Setting coordinates data file...")
    if user.use_rra_kinematics:
        kinfile = os.path.join(fpath, user.rracode, trial + "_RRA_" + str(user.rraiter) + "_Kinematics_q.sto")
        filtfreq = -1
    else:
        kinfile = os.path.join(fpath, user.ikcode, trial + "_ik.mot")    
        filtfreq = user.kinematics_filter_cutoff
    tool.setDesiredKinematicsFileName(kinfile)
    tool.setLowpassCutoffFrequency(filtfreq)


    
    # ******************************
    # RUN TOOL
    
    print("Running the CMCTool...\n---> (this may take a while)")
 
    # save the settings in a setup file
    customsetupfile = os.path.join(fpath, trial + "_Setup_CMC.xml")
    tool.printToXML(customsetupfile)   
 
    # run the tool (need to load the setup again into a new CMCTool and set the
    # use_fast_target flag)
    try:
        tool2 = opensim.CMCTool(customsetupfile)
        tool2.run()
        print("Done.")
    except:
        print("---> ERROR: CMC failed. Skipping CMC for %s." % trial)
    finally:
        shutil.copyfile(user.logfile, os.path.join(fpath, user.triallogfolder, "out_CMC.log"))
        print("------------------------------------------------\n")
    
    # ******************************
        
    return None    
    
    

'''
run_opensim_jr(osimkey, user):
    Set up and run the Tool using the API. A generic XML setup file is
    initially loaded, and then modified using the API. The Tool is then run
    via the API. Results are printed to text files in the remote folder.
'''
def run_opensim_jr(osimkey, user):
    
    # trial folder, model and trial
    fpath = osimkey.outpath
    modelfile = osimkey.model
    trial = osimkey.trial

    # clear log file
    open(user.logfile, "w").close()
    
    print("\nPerforming JR on trial: %s" % trial)
    print("------------------------------------------------")
    
    # create an generic AnalyzeTool from the template setup file (set the flag
    # aLoadModelAndInput = false as we are only configuring the setup file at 
    # this stage)
    print("Create new AnalyzeTool...")
    refsetuppath = user.refsetuppath
    refsetupfile = user.refsetupso
    tool = opensim.AnalyzeTool(os.path.join(refsetuppath, refsetupfile), False)
    tool.setName(trial) 
  
    # set the initial and final times (limit to between first and last event)
    t0 = float(osimkey.events["time"][0])
    t1 = float(osimkey.events["time"][osimkey.events["opensim_last_event_idx"]])
    print("Setting the time window: %0.3f sec --> %0.3f sec..." % (t0, t1))
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)
    
    # set output directory
    print("Setting output file name...")
    stofilepath = os.path.join(fpath, user.jrcode)
    if not os.path.isdir(stofilepath): os.makedirs(stofilepath)
    tool.setResultsDir(stofilepath)
    
    # create an external loads object from template, set it up, and
    # print to destination folder (This is not the best way to do this,
    # but Opensim Tools are designed to read directly from XML files,
    # so better to fully set up an external loads file, print it then
    # load it again into the Tool, than to create an ExternalLoads
    # object and connect it to the Model. This also ensures a copy of
    # the external loads file is available in the trial folder in case
    # a one-off analysis needs to be run in future.)
    print("Creating external loads XML file...")
    extloadsfile = os.path.join(fpath, trial + "_ExternalLoads.xml")
    if not os.path.isfile(extloadsfile):
        extloads = opensim.ExternalLoads(os.path.join(refsetuppath, user.additionalfilesfolder, user.refexternalloads), True)       
        extloads.setDataFileName(os.path.join(fpath, trial + "_grf.mot"))
        extloads.printToXML(extloadsfile)  

    # set the external loads file name
    tool.setExternalLoadsFileName(extloadsfile)
    
    
    # ******************************
    # CREATE JOINT REACTION ANALYSIS   
    
    # create an JR analysis
    print("Create new JR Analysis to the AnalysisSet...")
    jr = opensim.JointReaction()
    jr.setName("jr")
    jr.setStartTime(t0)
    jr.setEndTime(t1)
    
    # set the muscle forces 
    print("Setting muscle forces...")
    if user.jr_use_cmc_forces:
        forcefile = os.path.join(fpath, user.cmccode, trial + "_Actuation_force.sto")
    else:
        forcefile = os.path.join(fpath, user.socode, trial + "_so_force.sto")         
    jr.setForcesFileName(forcefile)
    
    # set the joints to be analysed
    print("Setting joints to be analysed...")
    joints = opensim.ArrayStr()
    onbodys = opensim.ArrayStr()
    inframes = opensim.ArrayStr()
    for j in user.jr_joints.keys():
        joints.append(j)
        onbodys.append(user.jr_joints[j][0])
        inframes.append(user.jr_joints[j][1])
    jr.setJointNames(joints)
    jr.setOnBody(onbodys)
    jr.setInFrame(inframes)
    
    # add the JointReaction analysis to the tool AnalysisSet
    print("Append new JR Analysis to the AnalysisSet...")
    analyses = tool.getAnalysisSet()
    analyses.insert(0, jr)


    # ******************************
    # PREPARE RESERVE ACTUATORS
    # (set the pelvis reserves to act at pelvis COM, or use default application
    # at the ground-pelvis origin located at pelvis geometric centre)
        
    # get the pelvis COM location from the model
    model = opensim.Model(os.path.join(fpath, modelfile))
    pelvis = opensim.Body.safeDownCast(model.findComponent("pelvis"))
    pelviscom = pelvis.getMassCenter() 
    
    # load reference actuator forceset, get the pelvis actuators and set force
    # application point to pelvis COM
    refrefreserveactuators = user.refreserveactuators
    residforceset = opensim.ForceSet(os.path.join(refsetuppath, user.additionalfilesfolder, refrefreserveactuators))
    for x in ["FX","FY","FZ"]:
        residforce = opensim.PointActuator.safeDownCast(residforceset.get(x))
        residforce.set_point(pelviscom)

    # write updated actuator set to file    
    forcesetfile = trial + "_Reserve_Actuators.xml"
    residforceset.printToXML(os.path.join(fpath, forcesetfile))
      
    # append reserve actuators
    # (note: need to create a new ArrayStr with one element, add the file name,
    # and then pass the ArrayStr to setForceSetFiles; however, the tool does
    # not like whitespace in the <force_set_files> XML property)
    print("Appending reserve actuators...")
    fsvec = opensim.ArrayStr()    
    fsvec.append(forcesetfile)
    tool.setForceSetFiles(fsvec)
    tool.setReplaceForceSet(False)   

    # # append default reference actuators, with pelvis residuals acting at 
    # # ground-pelvis origin located at pelvis geometric centre
    # refforcesetfile = user.refreserveactuators
    # residforceset = opensim.ForceSet(os.path.join(refsetuppath, refforcesetfile))
    # forcesetfile = trial + "_Reserve_Actuators.xml"
    # residforceset.printToXML(os.path.join(fpath, forcesetfile))
    # fsvec = opensim.ArrayStr()
    # fsvec.append(forcesetfile)
    # tool.setForceSetFiles(fsvec)
    # tool.setReplaceForceSet(False)     


    # ******************************
    # SET THE MODEL AND KINEMATICS
    # (use baseline or RRA adjusted model and kinematics)
           
    # desired model name
    print("Loading the model...")
    if user.jr_use_cmc_forces:
        if user.update_mass:
            actualmodelfile = trial + "_RRA_" + str(user.rraiter) + "_TorsoAdjusted_MassUpdated.osim"
        else:
            actualmodelfile = trial + "_RRA_" + str(user.rraiter) + "_TorsoAdjusted.osim"
    else:
        actualmodelfile = modelfile
    tool.setModelFilename(os.path.join(fpath, actualmodelfile))
        
    # set coordinates data file
    print("Setting coordinates data file...")
    if user.jr_use_cmc_forces:
        kinfile = os.path.join(fpath, user.cmccode, trial + "_Kinematics_q.sto")
        filtfreq = -1    
    else:            
        kinfile = os.path.join(fpath, user.ikcode, trial + "_ik.mot")
        filtfreq = 6.0
    tool.setCoordinatesFileName(kinfile)
    tool.setLowpassCutoffFrequency(filtfreq)
    
    
    # ******************************
    # RUN TOOL 
    
    print("Running the AnalysisTool (JR)...")

    # save the settings in a setup file
    customsetupfile = os.path.join(fpath, trial + "_Setup_JR.xml")
    tool.printToXML(customsetupfile)
    
    # run the tool (need to load the setup again into a new AnalyzeTool)
    try:
        tool2 = opensim.AnalyzeTool(customsetupfile)
        tool2.run()
        print("Done.")
    except:
        print("---> ERROR: JR failed. Skipping JR for %s." % trial)
    finally:
        shutil.copyfile(user.logfile, os.path.join(fpath, user.triallogfolder, "out_JR.log")) 
        print("------------------------------------------------\n")

    # ******************************
        
    return None  

    
    



'''
-----------------------------------
---- FUNCTIONS: MISCELLANEOUS -----
-----------------------------------
'''



'''
update_osim_fom(modelfullpath, scalefactor, refmodel):
    Update all OpenSim model FoM by a fixed scale factor, or using the scaling
    law described by Handsfield et al. 2013 (scalefactor["all"] = -1), and/or 
    apply a custom scale factor to selected muscles. To ignore scale factors,
    set scalefactor["all"] = 0.
    
    Usage:
        
        dict ("all", 0): do not scale (use model as-is)
        dict ("all", -1): apply Handsfield scaling law
        dict ("all", float): apply fixed scale factor to all muscles, this can
                be overwritten for specific muscles by adding additional dict
                items for each muscle as per below.
        dict (key, float): for each selected muscle (key) in the dict, apply
                custom scale factor (float). The muscle (key) may be a full
                muscle name (e.g. "vasint_r") or just the prefix (e.g.
                "vasint"). If the latter, then the scale factor is applied to
                all muscles found with that prefix.
            e.g. scalefactor = {}
                 scalefactor["all"] = 1.2
                 scalefactor["vasint"] = 1.5
                 scalefactor["semimem"] = 2.0        
'''
def update_osim_fom(model, scalefactor, refmodel):
    
    # load the model and get the muscles
    allmuscles = model.getMuscles()
    
    # scale by a fixed scale factor, or use Handsfield scaling law
    sf = 1.0
    if "all" in scalefactor:
        if scalefactor["all"] > 0:
            sf = scalefactor["all"]
            for m in range(allmuscles.getSize()):
                currmuscle = allmuscles.get(m)
                currmuscle.setMaxIsometricForce(sf * currmuscle.getMaxIsometricForce())
        elif scalefactor["all"] == -1:
            # Handsfield scaling law: TBD, requires refmodel
            return model
        elif scalefactor["all"] == 0:
            return model
                        
    # custom scale selected variables, overwrites muscles scaled by "all" key
    for sfname in scalefactor:
        if sfname.casefold() == "all": continue
        for m in range(allmuscles.getSize()):
            currmuscle = allmuscles.get(m)
            mname = currmuscle.getName()
            if mname.startswith(sfname):
                sfm = scalefactor[sfname] / sf
                currmuscle.setMaxIsometricForce(sfm * currmuscle.getMaxIsometricForce())

        
    return model



'''
update_osim_lom(modelfullpath, scalefactor):
    Update all OpenSim model LoM by a fixed scale factor, and/or apply a
    custom scale factor to selected muscles. To ignore scale factors, set
    scalefactor["all"] = 0.
    
    Usage:
        
        dict ("all", 0): do not scale (use model as-is)
        dict ("all", float): apply fixed scale factor to all muscles, this can
                be overwritten for specific muscles by adding additional dict
                items for each muscle as per below.
        dict (key, float): for each selected muscle (key) in the dict, apply
                custom scale factor (float). The muscle (key) may be a full
                muscle name (e.g. "vasint_r") or just the prefix (e.g.
                "vasint"). If the latter, then the scale factor is applied to
                all muscles found with that prefix.
            e.g. scalefactor = {}
                 scalefactor["all"] = 1.2
                 scalefactor["vasint"] = 1.5
                 scalefactor["semimem"] = 2.0
'''
def update_osim_lom(model, scalefactor):
    
    # load the model and get the muscles
    allmuscles = model.getMuscles()
    
    # scale by a fixed scale factor
    if "all" in scalefactor:
        if scalefactor["all"] > 0:
            sf = scalefactor["all"]
            for m in range(allmuscles.getSize()):
                currmuscle = allmuscles.get(m)
                currmuscle.setOptimalFiberLength(sf * currmuscle.getOptimalFibreLength())
        elif scalefactor["all"] == 0:
            return model
                        
    # custom scale selected variables, overwrites muscles scaled by "all" key
    for sfname in scalefactor:
        if sfname.casefold() == "all": continue
        for m in range(allmuscles.getSize()):
            currmuscle = allmuscles.get(m)
            mname = currmuscle.getName()
            if mname.startswith(sfname):
                sfm = scalefactor[sfname] / sf
                currmuscle.setOptimalFiberLength(sfm * currmuscle.getOptimalFiberLength())

        
    return model
            


'''
update_osim_lom_const_lmt(modelfullpath, scalefactor):
    Update all OpenSim model LoM by a fixed scale factor, and/or apply a
    custom scale factor to selected muscles. To ignore scale factors, set
    scalefactor["all"] = 0. Then adjust the LsT so that the total MTU rest
    length (LMT = LoM + LsT) remains constant.
    
    Usage:
        
        dict ("all", 0): do not scale (use model as-is)
        dict ("all", float): apply fixed scale factor to all muscles, this can
                be overwritten for specific muscles by adding additional dict
                items for each muscle as per below.
        dict (key, float): for each selected muscle (key) in the dict, apply
                custom scale factor (float). The muscle (key) may be a full
                muscle name (e.g. "vasint_r") or just the prefix (e.g.
                "vasint"). If the latter, then the scale factor is applied to
                all muscles found with that prefix.
            e.g. scalefactor = {}
                 scalefactor["all"] = 1.2
                 scalefactor["vasint"] = 1.5
                 scalefactor["semimem"] = 2.0
'''
def update_osim_lom_const_lmt(model, scalefactor):
    
    # load the model and get the muscles
    allmuscles = model.getMuscles()
    
    # scale all by a fixed scale factor, but also adjust the tendon slack
    # length to keep the total MTU rest length constant
    sf = 1.0
    if "all" in scalefactor:
        if scalefactor["all"] > 0:
            sf = scalefactor["all"]
            for m in range(allmuscles.getSize()):
                
                # current MTU properties
                currmuscle = allmuscles.get(m)                                
                lom0 = currmuscle.get_optimal_fiber_length()
                lst0 = currmuscle.get_tendon_slack_length()
                lmt = lom0 + lst0
                
                # new MTU properties
                lom1 = lom0 * sf
                lst1 = lmt - lom1
                currmuscle.set_optimal_fiber_length(lom1)
                currmuscle.set_tendon_slack_length(lst1)
                
        elif scalefactor["all"] == 0:
            return model

    # custom scale selected variables, adjust tendon slack length to ensure
    # total MTU rest length is constant     
    for sfname in scalefactor:
        if sfname.casefold() == "all": continue
        for m in range(allmuscles.getSize()):
            currmuscle = allmuscles.get(m)
            mname = currmuscle.getName()
            if mname.startswith(sfname):
                
                # current MTU properties
                lom0 = currmuscle.get_optimal_fiber_length()
                lst0 = currmuscle.get_tendon_slack_length()
                lmt = lom0 + lst0
                
                # new MTU properties
                sfm = scalefactor[sfname] / sf
                lom1 = lom0 * sfm
                lst1 = lmt - lom1
                currmuscle.set_tendon_slack_length(lom1)
                currmuscle.set_tendon_slack_length(lst1)
                    
    return model
            


'''
update_osim_lst(modelfullpath, scalefactor):
    Update all OpenSim model LsT by a fixed scale factor, and/or apply a
    custom scale factor to selected muscles. To ignore scale factors, set
    scalefactor["all"] = 0.
    
    Usage:
        
        dict ("all", 0): do not scale (use model as-is)
        dict ("all", float): apply fixed scale factor to all muscles, this can
                be overwritten for specific muscles by adding additional dict
                items for each muscle as per below.
        dict (key, float): for each selected muscle (key) in the dict, apply
                custom scale factor (float). The muscle (key) may be a full
                muscle name (e.g. "vasint_r") or just the prefix (e.g.
                "vasint"). If the latter, then the scale factor is applied to
                all muscles found with that prefix.
            e.g. scalefactor = {}
                 scalefactor["all"] = 1.2
                 scalefactor["vasint"] = 1.5
                 scalefactor["semimem"] = 2.0
'''
def update_osim_lst(model, scalefactor):
    
    # load the model and get the muscles
    allmuscles = model.getMuscles()
    
    # scale by a fixed scale factor
    if "all" in scalefactor:
        if scalefactor["all"] > 0:
            sf = scalefactor["all"]
            for m in range(allmuscles.getSize()):
                currmuscle = allmuscles.get(m)
                currmuscle.setTendonSlackLength(sf * currmuscle.getTendonSlackLength())
        elif scalefactor["all"] == 0:
            return model
                        
    # custom scale selected variables, overwrites muscles scaled by "all" key
    for sfname in scalefactor:
        if sfname.casefold() == "all": continue
        for m in range(allmuscles.getSize()):
            currmuscle = allmuscles.get(m)
            mname = currmuscle.getName()
            if mname.startswith(sfname):
                sfm = scalefactor[sfname] / sf
                currmuscle.setTendonSlackLength(sfm * currmuscle.getTendonSlackLength())

        
    return model



'''
write_ground_forces_mot_file(osimkey):
    Write .mot file for ground forces (external loads).
'''
def write_ground_forces_mot_file(osimkey):

    # output dataframe info
    ns = len(osimkey.forces["time"])
    nc = 19    

    # write headers
    fname = osimkey.trial + "_grf.mot"
    fpath = osimkey.outpath
    with open(os.path.join(fpath,fname),"w") as f:
        f.write("%s\n" % fname)
        f.write("version=1\n")
        f.write("nRows=%d\n" % ns)
        f.write("nColumns=%s\n" % nc)
        f.write("inDegrees=yes\n")
        f.write("endheader\n")

    # build data array
    datamat = np.zeros([ns, nc])
    datamat[:,0] = osimkey.forces["time"]
    datamat[:,1:4] = osimkey.forces["data"]["right"]["F"]
    datamat[:,4:7] = osimkey.forces["data"]["right"]["cop"]
    datamat[:,7:10] = osimkey.forces["data"]["left"]["F"]
    datamat[:,10:13] = osimkey.forces["data"]["left"]["cop"]
    datamat[:,13:16] = osimkey.forces["data"]["right"]["T"]
    datamat[:,16:19] = osimkey.forces["data"]["left"]["T"]    
        
    # convert to dataframe
    headers = ["time", "ground_force_vx", "ground_force_vy", "ground_force_vz", "ground_force_px", "ground_force_py", "ground_force_pz", "1_ground_force_vx", "1_ground_force_vy", "1_ground_force_vz", "1_ground_force_px", "1_ground_force_py", "1_ground_force_pz", "ground_torque_x", "ground_torque_y", "ground_torque_z", "1_ground_torque_x", "1_ground_torque_y", "1_ground_torque_z"]
    data = pd.DataFrame(datamat, columns=headers)
        
    # write table
    data.to_csv(os.path.join(fpath,fname), mode="a", sep="\t", header=True, index=False, float_format="%20.10f")
    
    return data



'''
write_marker_trajctory_trc_file(osimkey):
    Write .trc file for marker trajectories.
'''
def write_marker_trajctory_trc_file(osimkey):
    
    # output dataframe info
    ns = len(osimkey.markers["time"])
    nm = len(osimkey.markers) - 5
    nc = 2 + (nm * 3)
    rate = osimkey.markers["rate"]
    units = osimkey.markers["units"]

    # remove non-marker dict keys
    markernames0 = list(osimkey.markers.keys())
    markernames0.remove("rate")
    markernames0.remove("units")
    markernames0.remove("offset")
    markernames0.remove("frames")
    markernames0.remove("time")

    # build marker headers
    markernames = ""
    markernames = markernames.join(["Frame#\t","Time\t"] + list(map(lambda x: x + "\t\t\t", markernames0)))
    dirnums = list(range(1,len(markernames0) + 1))
    dirnames = ""
    dirnames = dirnames.join(["\t"] + list(map(lambda n: "\tX%d\tY%d\tZ%d" % (n, n, n), dirnums)))

    # write headers
    fname = osimkey.trial + "_markers.trc"
    fpath = osimkey.outpath
    with open(os.path.join(fpath,fname), "w") as f:
        f.write("PathFileType\t4\t(X/Y/Z)\t%s\n" % fname)
        f.write("DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n")
        f.write("%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n" % (rate, rate, ns, nm, "mm", rate, 1, ns))
        f.write("%s\n" % markernames)
        f.write("%s\n" % dirnames)
        
    # build data array
    datamat = np.zeros([ns, nc])
    datamat[:,0] = osimkey.markers["frames"]
    datamat[:,1] = osimkey.markers["time"]
    n = 2
    for mkr in markernames0:
        mkrdata = osimkey.markers[mkr]
        if units.casefold() == "m": mkrdata = mkrdata * 1000
        datamat[:,n:n+3] = mkrdata
        n = n + 3
    
    # convert to dataframe
    data = pd.DataFrame(datamat)
    data[0] = data[0].astype(int)
    
    # write table, no headers
    data.to_csv(os.path.join(fpath,fname), mode="a", sep="\t", header=False, index=False, float_format="%20.10f")
    
    return data



'''
perform_recommended_mass_change(rra_model, rra_log_file_fullpath):
    Extract the recommended mass changes from the RRA log file and update the
    model segment masses. The parameter rra_log_file_fullpath should contain 
    the full path to the relevant log file.
'''
def perform_recommended_mass_change(rra_model, rra_log_file_fullpath):
    
    # define the output format
    mstr = re.compile("\*  (\w+): orig mass = (\d+.\d+), new mass = (\d+.\d+)")
        
    # parse the log file and extract current and recommended masses
    new_mass_set = {}
    with open(rra_log_file_fullpath, "r") as f:
        for line in f:
            tokens = mstr.match(line)
            if tokens:
                new_mass_set[tokens[1]] = {}
                new_mass_set[tokens[1]]["old"] = float(tokens[2])
                new_mass_set[tokens[1]]["new"] = float(tokens[3])              
    
    # update segment masses
    bodyset = rra_model.getBodySet()
    for key in new_mass_set:
        body = bodyset.get(key)
        body.setMass(new_mass_set[key]["new"])
               
    return rra_model



'''
prescribe_kinematics(model, ikdata, coord_list):
    Prescribe the coordinates listed in coord_list using the data from 
    ikdata (which can be IK results or RRA results).
'''
def prescribe_kinematics(model, ikdata, coord_list, butter_order, cutoff):
    
    # get the time vector
    timevec = ikdata.get("time").to_list()
    
    # calculate the sample rate (this can be different to the original C3D
    # sample rate if RRA kinematics are used instead of IK)
    sample_rate = 1 / (timevec[1] - timevec[0])
    
    # create prescribed joint kinematics object and add to model
    coordset = model.getCoordinateSet()
    for coordname in coord_list:
        
        # x, y data series
        xvals = timevec
        yvals = ikdata.get(coordname).to_list()
        
        # if not pelvis transation, convert to radians
        if not coordname.startswith("pelvis_t"):
            yvals = np.radians(yvals).tolist()
        
        # filter data
        yvals = filter_timeseries(yvals, sample_rate, butter_order, cutoff)
        
        # create a SimmSpline for the timeseries
        prescribed_spline = opensim.SimmSpline()
        for n in range(1,len(timevec)):
            prescribed_spline.addPoint(xvals[n], yvals[n])
            
        # apply SimmSpline to coordinate
        coordinate = coordset.get(coordname)
        coordinate.set_prescribed(True)
        coordinate.setPrescribedFunction(prescribed_spline)
        
    return model
   
         
        
'''
filter_timeseries(data_raw, sample_rate, butter_order, cutoff):
    Filter timeseries data. Raw data can be a list, or an array with rows
    representing time steps and columns as variables.
    
    Duplicates the same function in c3dextract.py.
''' 
def filter_timeseries(data_raw, sample_rate, butter_order, cutoff):
      
      # filter design
      Wn = sample_rate / 2
      normalised_cutoff = cutoff / Wn
      b, a = signal.butter(butter_order, normalised_cutoff, "lowpass")
      
      # apply filter
      data_filtered = signal.filtfilt(b, a, data_raw, axis = 0)

      return data_filtered      
