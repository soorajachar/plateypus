#!/usr/bin/env python3 
import json,pickle,os
import numpy as np
import pandas as pd
import tkinter as tk
from plateypus.dataprocessing import initialDataProcessing as idp
from plateypus.dataprocessing import cytokineDataProcessing as cydp
from plateypus.dataprocessing import cellDataProcessing as cdp
from plateypus.dataprocessing import proliferationDataProcessing as pdp
from plateypus.dataprocessing import killingDataProcessing as kdp
from plateypus.dataprocessing import singleCellDataProcessing as scdp
#import automatedCBAProcessingGUI as autoCBA


secondPath = '../../outputData'

concUnit = 1e9
unitPrefixDictionary = {1e12:'pM',1e9:'nM',1e6:'uM',1e3:'mM',1e0:'M'}
concUnitPrefix = unitPrefixDictionary[concUnit]

class DataProcessingStartPage(tk.Frame):
    def __init__(self, master,folderName,expNum,ex_data,bPage):
        tk.Frame.__init__(self, master)
        
        #os.chdir(master.homedirectory+'/'+folderName)
        backPage = bPage

        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l2 = tk.Label(mainWindow, text="""Datatype: """,padx = 20).grid(row=0,column=0,sticky=tk.W)
        l2a = tk.Label(mainWindow, text="Cytokine:",padx = 20).grid(row=1,column=0,sticky=tk.W)
        l2b = tk.Label(mainWindow,text="Cell:",padx = 20).grid(row=2,column=0,sticky=tk.W)
        l2c = tk.Label(mainWindow,text="Proliferation:",padx = 20).grid(row=3,column=0,sticky=tk.W)
        l2e = tk.Label(mainWindow,text="Killing:",padx = 20).grid(row=4,column=0,sticky=tk.W)
        l2d = tk.Label(mainWindow,text="Single Cell:",padx = 20).grid(row=5,column=0,sticky=tk.W)
        
        def createDataFrame(dataType):
            dataProcessingMaster(folderName,expNum,dataType,ex_data,v3.get())
        
        l3 = tk.Label(mainWindow, text="""Action: """).grid(row=0,column=1,sticky=tk.W)
        
        cytCalibrationParametersButton = tk.Button(mainWindow,text='Enter CBA bead calibration parameters',command=lambda: master.switch_frame(cydp.CalibrationParameterPage,folderName,expNum,ex_data,DataProcessingStartPage,bPage))
        cytCalibrationParametersButton.grid(row=1,column=1,sticky=tk.W)
        #cytCalibrationParametersButton = tk.Button(mainWindow,text='Create CBA gates',command=lambda: master.switch_frame(autoCBA.AutomateCBAStartPage,folderName,expNum,ex_data,DataProcessingStartPage,bPage))
        #cytCalibrationParametersButton.grid(row=1,column=2,sticky=tk.W)
        cytDfButton = tk.Button(mainWindow,text='Create dataframe',command=lambda: createDataFrame('cyt'))
        cytDfButton.grid(row=1,column=2,sticky=tk.W)

        cellDfButton = tk.Button(mainWindow,text='Create dataframe',command=lambda: createDataFrame('cell'))
        cellDfButton.grid(row=2,column=2,sticky=tk.W)
        #cellAbPanelButton = tk.Button(mainWindow,text='Edit antibody panel',command=lambda: master.switch_frame(cdp.MarkerNumberPage,folderName,expNum,ex_data,DataProcessingStartPage,bPage))
        #cellAbPanelButton.grid(row=2,column=2,sticky=tk.W)
        
        prolifGenerationGatesButton = tk.Button(mainWindow,text='Edit generation gates',command=lambda: master.switch_frame(pdp.ProliferationSelectionPage,folderName,expNum,ex_data,DataProcessingStartPage,backPage))
        prolifGenerationGatesButton.grid(row=3,column=1,sticky=tk.W)
        prolifDfButton = tk.Button(mainWindow,text='Create dataframe',command=lambda: createDataFrame('prolif'))
        prolifDfButton.grid(row=3,column=2,sticky=tk.W)
        
        killingDfButton = tk.Button(mainWindow,text='Create dataframe',command=lambda: createDataFrame('killing'))
        killingDfButton.grid(row=4,column=2,sticky=tk.W)

        #completeSingleCellDfButton = tk.Button(mainWindow,text='Create complete dataframes')
        #completeSingleCellDfButton.grid(row=4,column=2,sticky=tk.W)
        singleCellDfButton = tk.Button(mainWindow,text='Create dataframe',command=lambda: createDataFrame('singlecell'))
        singleCellDfButton.grid(row=5,column=2,sticky=tk.W)
        cbWindow = tk.Frame(mainWindow)
        cbWindow.grid(row=5,column=1,sticky=tk.W)
        l3 = tk.Label(cbWindow,text='Use empty wells?').grid(row=0,column=0,sticky=tk.W)
        v3 = tk.BooleanVar(value=False)
        cb = tk.Checkbutton(cbWindow, variable=v3)
        cb.grid(row=0,column=1,sticky=tk.W)

        for i,button in enumerate([cytDfButton,cellDfButton,prolifDfButton,prolifGenerationGatesButton,singleCellDfButton,killingDfButton]):
            if i == 0:
                requiredFiles = ['CBAcalibrationParameters-'+folderName+'.json']
            elif i == 2:
                requiredFiles = ['singleCellDataFrame-proliferation-'+folderName+'.pkl']
            #Go to separate page where you select how many markers are proliferation, then each marker that serves as proliferation
            elif i == 3:
                requiredFiles = ['initialSingleCellDf-channel-'+folderName+'.h5']
                #requiredFiles = ['logicleProliferationDf.pkl','rawProliferationDf.pkl']
            else:
                requiredFiles = []
            for requiredFile in requiredFiles:
                if requiredFile not in os.listdir('misc')+os.listdir('inputData/bulkCSVFiles')+os.listdir('outputData/pickleFiles'):
                    button.config(state=tk.DISABLED)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        #tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

def dataProcessingMaster(folderName,expNum,dataType,ex_data,useBlankWells):
    if dataType == 'singlecell' or dataType == 'prolif':
        parameterExtension = 'cell'
    else:
        parameterExtension = dataType
    experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'-'+parameterExtension+'.json','r'))
    if experimentParameters['format'] == 'plate':
        experimentFormat = 'plate'
        experimentLevelLayoutDict = pickle.load(open('misc/layoutDict-'+folderName+'-'+parameterExtension+'.pkl','rb'))
    else:
        experimentFormat = 'tube'
        experimentLevelLayoutDict = pickle.load(open('misc/tubeLayout-'+folderName+'-'+parameterExtension+'.pkl','rb'))
    #experimentLevelLayoutDict = idp.tilePlateLayouts(experimentParameters,levelLayouts)
    if(dataType == 'cyt'):
        calibrationParameters = json.load(open('misc/CBAcalibrationParameters-'+folderName+'.json','r'))
        numberOfCalibrationSamples = calibrationParameters['Number']
        initialStandardVolume = calibrationParameters['Volume']
        cydp.calibrateExperiment(folderName,secondPath,concUnit,concUnitPrefix,numberOfCalibrationSamples,initialStandardVolume)
        basecytdf = idp.createBaseDataFrame(experimentParameters,folderName,expNum,dataType,experimentLevelLayoutDict)
        cytdf = cydp.createCytokineDataFrame(folderName,basecytdf,concUnitPrefix)
        idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,cytdf,ex_data) 
    elif(dataType == 'cell'):
        celldf = idp.createBaseDataFrame(experimentParameters,folderName,expNum,dataType,experimentLevelLayoutDict)
        idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,celldf,ex_data) 
    elif(dataType == 'prolif'):
        prolifdf = pdp.generateBulkProliferationStatistics(folderName,expNum)
        idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,prolifdf,ex_data) 
    elif(dataType == 'killing'):
        killingdf = kdp.generateBulkKillingStatistics(experimentParameters,folderName,expNum,dataType,experimentLevelLayoutDict)
        idp.saveFinalDataFrames(folderName,secondPath,expNum,dataType,killingdf,ex_data) 
        
        #Re-add parsed time values back into experiment parameter file
        #levelLabelDict = experimentParameters['levelLabelDict']
        #levelLabelDict['Time'] = list(killingDf.columns) 
        #experimentParameters['levelLabelDict'] = levelLabelDict
        #with open('misc/experimentParameters-'+folderName+'-'+dataType+'.json', 'w') as fp:
        #    json.dump(experimentParameters, fp)

    elif(dataType == 'singlecell'):
        dataType = 'singlecell'
        if experimentFormat == 'plate':
            scdp.createPlateSingleCellDataFrame(folderName,experimentParameters,experimentLevelLayoutDict,useBlankWells)
        else:
            scdf = scdp.createTubeSingleCellDataFrame(folderName,experimentParameters,experimentLevelLayoutDict)
    niceDatatypeNameDict = {'cyt':'Cytokine','cell':'Bulk cell','prolif':'Proliferation','singlecell':'Single cell','killing':'Killing'}
    tk.messagebox.showinfo("Dataframe Created", niceDatatypeNameDict[dataType]+" dataframe created!")
