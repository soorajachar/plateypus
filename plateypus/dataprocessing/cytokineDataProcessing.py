#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,glob
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import colorcet as cc
import numpy as np
import pandas as pd
import seaborn as sns
import tkinter as tk
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from itertools import groupby
from tkinter.ttk import Combobox
from tkinter.messagebox import askyesno,showinfo
import os
import warnings,shutil
warnings.filterwarnings("ignore")
if os.name == 'nt':
    dirSep = '\\'
else:
    dirSep = '/'
from plateypus.dataprocessing.miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV, setMaxWidth
idx = pd.IndexSlice

class CytokineParsingPage(tk.Frame):
    def __init__(self, master,completeCytokineMWDf,unparsedCytokines,folderName,expNum,ex_data,shp,bPage):
        tk.Frame.__init__(self, master)
        
        secondaryhomepage = shp
        
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text='Several cytokine names in raw data were not found in loaded CBA kit.\nPlease choose their correct names below:').pack()
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10) 
        
        tk.Label(mainWindow, text='Unparsed Name:').grid(row=0, column=0)
        tk.Label(mainWindow, text='Kit Name:').grid(row=0, column=1)
        
        kitCytokineNames = list(pd.unique(completeCytokineMWDf.reset_index()['Cytokine']))
        kitCytokineNameVarList = []
        for i,unparsedCytokine in enumerate(unparsedCytokines):
            tk.Label(mainWindow, text=unparsedCytokine).grid(row=i+1, column=0)
            kitCytokineNameVar = tk.StringVar()
            kitCytokineNameDropdown = tk.OptionMenu(mainWindow,kitCytokineNameVar,*kitCytokineNames)
            kitCytokineNameDropdown.grid(row=i+1, column=1)
            kitCytokineNameVarList.append(kitCytokineNameVar)
         
        def collectInputs():
            cytokineRenamingDict = {x:y.get() for x,y in zip(unparsedCytokines,kitCytokineNameVarList)}
            with open('misc'+dirSep+'cytRenamingDict-'+folderName+'.json','w') as f:
                json.dump(cytokineRenamingDict,f)
            master.switch_frame(secondaryhomepage,folderName,expNum,ex_data,bPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class CalibrationParameterPage(tk.Frame):
    def __init__(self, master,folderName,expNum,ex_data,shp,bPage):
        tk.Frame.__init__(self, master)
        
        secondaryhomepage = shp
        
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10) 
        
        l1 = tk.Label(mainWindow,text='Number of CBA calibration samples: ').grid(row=0,column=0)
        t1 = tk.Entry(mainWindow)
        t1.grid(row=1,column=0,sticky=tk.W)

        l2 = tk.Label(mainWindow,text='Volume of initial CBA calibration solution (mL): ').grid(row=0,column=1)
        t2 = tk.Entry(mainWindow)
        t2.grid(row=1,column=1,sticky=tk.W)
        
        l4 = tk.Label(mainWindow,text='Serial Dilution Factor: ').grid(row=0,column=2)
        t4 = tk.Entry(mainWindow)
        t4.grid(row=1,column=2,sticky=tk.W)
        
        l3 = tk.Label(mainWindow, text='Species of CBA samples: ').grid(row=0, column=3)

        kitDf = pd.read_csv(master.homedirectory+'misc'+dirSep+'kitDf.csv').set_index(['Species','Kit'])
        self.kitDict = {x:kitDf.query("Species == @x").index.unique('Kit').tolist() for x in kitDf.index.unique('Species').tolist()}
        
        speciesVar = tk.StringVar(self)
        kitVar = tk.StringVar(self)
        
        kitMenu = tk.OptionMenu(mainWindow, kitVar, '')
        kitMenu.grid(row=1, column=4)
        
        def update_options(*args):
            kits = self.kitDict[speciesVar.get()]
            kitVar.set(kits[0])
            
            menu = kitMenu['menu']
            menu.delete(0, 'end')

            for kit in kits:
                menu.add_command(label=kit, command=lambda k=kit: kitVar.set(k))
            
        speciesVar.trace('w', update_options)
        speciesVar.set('Mouse')
        
        speciesMenu = tk.OptionMenu(mainWindow, speciesVar, *self.kitDict.keys())
        speciesMenu.grid(row=1, column=3)

        l4 = tk.Label(mainWindow, text='Cytokine kit used:').grid(row=0, column=4)
        
        tk.Button(mainWindow, text="Create new kit",command=lambda: master.switch_frame(KitCreationPage,folderName,expNum,ex_data,secondaryhomepage,bPage)).grid(row=3,column=4)
         
        def collectInputs():
            numCalibrationSamples = int(t1.get())
            initialStandardVolume = float(t2.get())
            kit = kitVar.get()
            sdf = float(t4.get())
            calibrationParameterDict = {'Number': numCalibrationSamples,'Volume': initialStandardVolume, 'Kit':[kit], 'SerialDilutionFactor':sdf, 'Species':speciesVar.get()}
            with open('misc'+dirSep+'CBAcalibrationParameters-'+folderName+'.json','w') as f:
                json.dump(calibrationParameterDict,f)
            master.switch_frame(secondaryhomepage,folderName,expNum,ex_data,bPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class KitCreationPage(tk.Frame):
    def __init__(self, master,folderName,expNum,ex_data,shp,bPage):
        tk.Frame.__init__(self, master)

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10) 
        
        tk.Label(mainWindow, text='Species: ').grid(row=0, column=0)
        tk.Label(mainWindow, text='Cytokine: ').grid(row=1, column=0)
        tk.Label(mainWindow, text='MW (Daltons): ').grid(row=2,column=0)
        tk.Label(mainWindow, text='Mass (pg): ').grid(row=3,column=0)

        #cytDf = pd.read_pickle(master.homedirectory+'misc'+dirSep+'cytokineMWDf.pkl').set_index('Cytokine',append=True)
        cytDf = pd.read_csv(master.homedirectory+'misc'+dirSep+'cytokineMWDf.csv').set_index(['Species','Cytokine'])
        self.cytokineDict = {x:cytDf.query("Species == @x").index.unique('Cytokine').tolist() for x in cytDf.index.unique('Species').tolist()}
        self.mwDict = {x+'--'+y:z for x,y,z in zip(list(cytDf.index.get_level_values('Species')),list(cytDf.index.get_level_values('Cytokine')),cytDf['MW'])}
        
        self.cb1 = Combobox(mainWindow, state="readonly", values=list(self.cytokineDict.keys()))
        self.cb1.current(0)
        self.cb1.grid(row=0,column=1)
        
        kitButtonWindow = tk.Frame(mainWindow)
        kitButtonWindow.grid(row=5,column=0,columnspan=2,pady=(20,0))
        
        #kitDf = pd.read_pickle(master.homedirectory+'misc'+dirSep+'kitDf.pkl')
        kitDf = pd.read_csv(master.homedirectory+'misc'+dirSep+'kitDf.csv').set_index(['Species','Kit'])
        species = self.cb1.get()
        self.kitNameCombo = Combobox(kitButtonWindow,values=kitDf.query("Species == @species").index.unique('Kit').tolist())
        tk.Label(kitButtonWindow,text='Kit: ').grid(row=0,column=0)
        self.kitNameCombo.grid(row=0,column=1)

        def populate_slave():
            self.cb2['values'] = self.cytokineDict[self.cb1.get()]
            self.cb2.current(0)
            populate_slave2()
            
            species = self.cb1.get()
            #kitDf = pd.read_pickle(master.homedirectory+'misc'+dirSep+'kitDf.pkl')
            kitDf = pd.read_csv(master.homedirectory+'misc'+dirSep+'kitDf.csv').set_index(['Species','Kit'])
            self.kitNameCombo['values'] = kitDf.query("Species == @species").index.unique('Kit').tolist() 
            
            return True #must return True or validatecommand will quit
        
        self.cb1.config(validate="focus", validatecommand=populate_slave)

        self.cb2 = Combobox(mainWindow)
        self.cb2.grid(column=1, row=1)
         
        self.mwEntry = tk.Entry(mainWindow,width=10)
        self.mwEntry.grid(row=2,column=1)
        
        self.massEntry = tk.Entry(mainWindow,width=5)
        self.massEntry.grid(row=3,column=1)
        self.massEntry.insert(0,str(10000))

        def populate_slave2():
            if self.cb1.get()+'--'+self.cb2.get() in self.mwDict:
                self.mwEntry.delete(0,tk.END)
                mw = self.mwDict[self.cb1.get()+'--'+self.cb2.get()]
                self.mwEntry.insert(0,str(mw))
            return True #must return True or validatecommand will quit

        self.cb2.config(validate="focus", validatecommand=populate_slave2)
        populate_slave()
        populate_slave2()
        
        #self.currentKitSpeciesList = []
        self.currentKitCytokineList = []
        self.currentKitMWList = []
        self.currentKitMassList = []
        def add_cytokine():
            #self.currentKitSpeciesList.append(self.cb1.get())
            if self.cb2.get() not in self.currentKitCytokineList:
                self.currentKitCytokineList.append(self.cb2.get())
                self.currentKitMWList.append(int(self.mwEntry.get()))
                self.currentKitMassList.append(int(self.massEntry.get()))
            else:
                cytIndex = self.currentKitCytokineList.index(self.cb2.get())
                self.currentKitCytokineList[cytIndex] = self.cb2.get()
                self.currentKitMWList[cytIndex] = int(self.mwEntry.get())
                self.currentKitMassList[cytIndex] = int(self.massEntry.get())

            print(self.cb1.get()+' '+self.cb2.get()+' with MW '+str(self.mwEntry.get())+' added\nCurrent cytokine list: '+','.join(self.currentKitCytokineList)+'\n')
        
        def remove_cytokine():
            #removedSpecies = self.currentKitSpeciesList.pop()
            removedSpecies = self.cb1.get() 
            removedCyt = self.currentKitCytokineList.pop()
            removedMW = self.currentMWList.pop()
            removedConc = self.currentConcList.pop()
            print(removedSpecies+' '+removedCyt+' with MW '+str(removedMW)+' removed\nCurrent cytokine list: '+','.join(self.currentKitCytokineList)+'\n')

        cytokineButtonWindow = tk.Frame(mainWindow)
        cytokineButtonWindow.grid(row=4,column=0,columnspan=2)
        
        tk.Button(cytokineButtonWindow, text="Add cytokine",command=lambda: add_cytokine()).grid(row=0,column=0)
        tk.Button(cytokineButtonWindow, text="Remove cytokine",command=lambda: remove_cytokine()).grid(row=0,column=1)

        kitButtonWindow2 = tk.Frame(mainWindow)
        kitButtonWindow2.grid(row=6,column=0,columnspan=2)
        def add_kit():
            #oldKitDf = pd.read_pickle(master.homedirectory+'misc'+dirSep+'kitDf.pkl')
            oldKitDf = pd.read_csv(master.homedirectory+'misc'+dirSep+'kitDf.csv').set_index(['Species','Kit'])
            listLength = len(self.currentKitCytokineList)
            newKitDf = pd.DataFrame({'Species':[self.cb1.get()]*listLength,'Kit':[self.kitNameCombo.get()]*listLength,'Cytokine':self.currentKitCytokineList,'MW':self.currentKitMWList,'Mass':self.currentKitMassList}).set_index(oldKitDf.index.names)
            kitDf = pd.concat([oldKitDf,newKitDf])
            tempDfList = []
            for s in sorted(kitDf.index.unique('Species').tolist()):
                sDf = kitDf.query("Species == @s")
                tempDf = pd.concat([sDf.query("Kit == @x") for x in sorted(sDf.index.unique('Kit').tolist())])
                tempDfList.append(tempDf)
            kitDf = pd.concat(tempDfList)
            additionalMessage = '\nThis kit contains these cytokines: \n'+','.join(self.currentKitCytokineList)
            answer = askyesno(title='Confirmation',message='Are you sure you want to add '+self.cb1.get()+' '+self.kitNameCombo.get()+' kit?'+additionalMessage)
            if answer:
                #kitDf.drop_duplicates().to_pickle(master.homedirectory+'misc'+dirSep+'kitDf.pkl')
                temp = kitDf.reset_index().drop_duplicates()
                dList = []
                for species in pd.unique(temp['Species']):
                    sDf  = temp.query("Species == @species")
                    for kit in pd.unique(sDf['Kit']):
                        kDf = sDf.query("Kit == @kit")
                        dList.append(kDf)
                temp = pd.concat(dList)
                temp.to_csv(master.homedirectory+'misc'+dirSep+'kitDf.csv',index=None)
                self.currentKitCytokineList = []
                self.currentKitMWList = []
                populate_slave()
                showinfo(title='Kit added',message=self.cb1.get()+' '+self.kitNameCombo.get()+' kit added!')

        def remove_kit():
            spec = self.cb1.get()
            k = self.kitNameCombo.get() 
            answer = askyesno(title='Confirmation', message='Are you sure you want to remove '+spec+' '+k+' kit?')
            if answer:
                #kitDf = pd.read_pickle(master.homedirectory+'misc'+dirSep+'kitDf.pkl')
                kitDf = pd.read_csv(master.homedirectory+'misc'+dirSep+'kitDf.csv').set_index(['Species','Kit'])
                #kitDf.query("Kit != @k").to_pickle(master.homedirectory+dirSep+'misc/kitDf.pkl')
                kitDf.query("Kit != @k").reset_index().to_csv(master.homedirectory+dirSep+'misc/kitDf.csv',index=None)
                populate_slave()
                showinfo(title='Kit removed',message=spec+' '+k+' kit removed!')

        tk.Button(kitButtonWindow2, text="Add kit",command=lambda: add_kit()).grid(row=0,column=0)
        tk.Button(kitButtonWindow2, text="Remove kit",command=lambda: remove_kit()).grid(row=0,column=1)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
         
        def collectInputs():
            answer = askyesno(title='Confirmation', message='Are you sure you want to exit cytokine kit creation?')
            if answer:
                master.switch_frame(shp,folderName,expNum,ex_data,bPage)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

def parseCytokineCSVHeaders(columns,cytokineMWDf,cytRenamingDict):
    #,Beads/IFNg | Geometric Mean (YG586-A),Beads/IL-2 | Geometric Mean (YG586-A),Beads/IL-4 | Geometric Mean (YG586-A),Beads/IL-6 | Geometric Mean (YG586-A),Beads/IL-10 | Geometric Mean (YG586-A),Beads/IL-17A | Geometric Mean (YG586-A),Beads/TNFa | Geometric Mean (YG586-A),
    newMultiIndexList = []
    for column in columns[1:-1]:
        populationNameVsStatisticSplit = column.split(' | ')
        if 'Barcode' not in populationNameVsStatisticSplit[0]:
            cytokine = populationNameVsStatisticSplit[0].split('/')[-1]
        else:
            cytokine = populationNameVsStatisticSplit[0].split('/')[-2]
        if cytokine in cytRenamingDict:
            cytokine = cytRenamingDict[cytokine]
        newMultiIndexList.append([cytokine])
    
    unparsedCytokines = []
    if type(cytokineMWDf) != list:
        kitCytokines = list(pd.unique(cytokineMWDf.reset_index()['Cytokine']))
        for cyt in newMultiIndexList:
            if cyt[0] not in kitCytokines:
                unparsedCytokines.append(cyt[0])

    return newMultiIndexList,unparsedCytokines

def performCommaCheck(fileName):
    with open('inputData'+dirSep+'bulkCSVFiles'+dirSep+fileName, 'r') as istr:
        with open('inputData'+dirSep+'bulkCSVFiles'+dirSep+'temp-'+fileName, 'w') as ostr:
            for line in istr:
                if line[-1] != ',':
                    line = line.rstrip('\n') + ','
                    print(line, file=ostr)
                else:
                    line = line.rstrip('\n')
                    print(line, file=ostr)
    os.remove('inputData'+dirSep+'bulkCSVFiles'+dirSep+fileName)
    os.rename('inputData'+dirSep+'bulkCSVFiles'+dirSep+'temp-'+fileName,'inputData'+dirSep+'bulkCSVFiles'+dirSep+fileName)

#Create Calibration Curves, obtain LODs (limits of detection) of experiment
def calibrateExperiment(folderName,secondPath,concUnit,concUnitPrefix,numberOfCalibrationSamples,initialStandardVolume,serialDilutionFactor,completeCytokineMWDf,cytRenamingDict={}):
    #Get cytokine calibration curve data
    tempExperimentParameters = {'overallPlateDimensions':[8,12]}
    calibrationFileNames = glob.glob('inputData'+dirSep+'bulkCSVFiles'+dirSep+'Calibration*')
    calibrationNames = []
    kitNames = []
    for calibrationFileName in calibrationFileNames:
        #performCommaCheck(calibrationFileName.split('/')[-1])
        newName = calibrationFileName.split('.')[0].split('_')[0].split(dirSep)[-1]
        kitNames.append(newName)
    sortedData,sortedFiles = cleanUpFlowjoCSV(kitNames,folderName,'cyt',tempExperimentParameters)
    for i,newName in enumerate(kitNames):
        if '-' in newName:
            newName2 = newName.split('-')[1]
        else:
            newName2 = newName
        kitNames[i] = newName2
    rsquaredList = []
    concLODList = []
    fittingParametersList = []
    cbaStandardsMFIList = []
    cbaPlotPointsMFIList = []
    cbaStandardsConcentrationList = []
    cbaPlotPointsConcentrationList = []

    numberOfPlotPoints = 101
    xaxistitle = 'Concentration of Cytokine Standards Standards ('+concUnitPrefix+')' 
    yaxistitle = 'GeoMFI'

    # Combine cytokine MWs into a dict
    completeCytokineMWDict = {x:y for x,y in zip(completeCytokineMWDf.reset_index()['Cytokine'],completeCytokineMWDf.reset_index()['MW'])}
    completeCytokineMassDict = {x:y for x,y in zip(completeCytokineMWDf.reset_index()['Cytokine'],completeCytokineMWDf.reset_index()['Mass'])}
    
    #TODO: Use fuzzy here to match mispelled cytokine names
    allCytokinesHaveMWInDict = True
    unparsedCytokinesBool = False
    for calibration in sortedData:
        cytokines,unparsedCytokines = parseCytokineCSVHeaders(calibration.columns,completeCytokineMWDf,cytRenamingDict)
        if len(unparsedCytokines) == 0:
            cytokines = [x[0] for x in cytokines]
            for cytokine in cytokines:
                if cytokine not in completeCytokineMWDict:
                    allCytokinesHaveMWInDict = False
        else:
            unparsedCytokinesBool = True
            break
    
    if not unparsedCytokinesBool:
        for calibration in sortedData:
            data = np.array(calibration.values[:,1:-1],dtype=float)
            cytokines,unparsedCytokines = parseCytokineCSVHeaders(calibration.columns,completeCytokineMWDf,cytRenamingDict)
            cytokines = [x[0] for x in cytokines]
            fittingParameters = np.zeros((data.shape[1],4))
            concLOD = np.zeros((data.shape[1],4))
            
            #masses in pg
            masses = np.array([completeCytokineMassDict[cyt] for cyt in cytokines])
            #conc in g/L
            initialConc = (masses*1e-12) /((initialStandardVolume*1e-3))
            
            #Calibration samples are always diluted by a factor of serialdilutionFactor (so with 12 calibration samples, the last sample is (serialDilutionFactor^-11) the concentration of the first, which is pure standard (2^0)
            x = np.linspace(-numberOfCalibrationSamples+1,0,numberOfCalibrationSamples)
            cbaStandardsConcentrations = np.fliplr(np.multiply(np.tile(initialConc, (x.shape[0],1)).T,np.power(serialDilutionFactor,np.tile(x, (initialConc.shape[0],1)))))
            #More x values along the above concentration bounds are sampled to use to construct calibration curve. Plot points are extended slightly at high range to allow visualization of upper LOD not accessible with experimental dilution
            x = np.linspace(-numberOfCalibrationSamples+1,4,numberOfPlotPoints)
            cbaStandardsConcentrationsPlotPoints = np.fliplr(np.multiply(np.tile(initialConc, (x.shape[0],1)).T,np.power(serialDilutionFactor,np.tile(x, (initialConc.shape[0],1)))))
            
            cbaStandardsConcentrationMatrix = np.zeros([len(cytokines),cbaStandardsConcentrations.shape[1]]) 
            cbaStandardsConcentrationPlotPointsMatrix = np.zeros([len(cytokines),cbaStandardsConcentrationsPlotPoints.shape[1]])
            cbaStandardsMFIMatrix = np.zeros([len(cytokines),cbaStandardsConcentrations.shape[1]]) 
            cbaStandardsMFIPlotPointsMatrix = np.zeros([len(cytokines),cbaStandardsConcentrationsPlotPoints.shape[1]])
            color_list = sns.color_palette(sns.color_palette(),len(cytokines))
            rsquared_kit = []
            for i,cytokineList in enumerate(cytokines):
                #amplitude bounded from range/2 to range*2, EC50 bounded from minimum to maximum standard concentration tested, Hill coefficient bounded from 0 to 2, Background bounded from 0 to minimum GFI*2
                lowerCurveFitBounds = [(np.max(data[:,i])-np.min(data[:,i]))/2,np.min(cbaStandardsConcentrations),0,0]
                upperCurveFitBounds = [(np.max(data[:,i])-np.min(data[:,i]))*2, np.max(cbaStandardsConcentrations), 2,np.min(data[:,i])*2]
                #use scipy curve fit to determine best hill equation fit for data, searching within the bounds given above
                popt,pcov = curve_fit(Hill, cbaStandardsConcentrations[i,:],np.log10(data[:,i]),sigma=np.log10(data[:,i]),bounds=(lowerCurveFitBounds,upperCurveFitBounds))
                rsquared = round(r_squared(cbaStandardsConcentrations[i,:],np.log10(data[:,i]),Hill,popt),3)
                rsquared_kit.append(rsquared)
                for j in range(len(popt)):  
                    #Convert just ec50 value to desired units (nM,uM etc) if cytokine has a molar mass in dict
                    if j == 1 and allCytokinesHaveMWInDict:
                        fittingParameters[i,j] = np.multiply(popt[j],(concUnit/completeCytokineMWDict[cytokine]))
                    #other values in 4 parameter logistic equation are tied to intensity y-value, which doesn't change, or are the hill coefficient, which is completely separate, so parameters are kept the same
                    else:
                        fittingParameters[i,j] = popt[j]
                
                #Convert x values of experimental data points and curve fit points to desired units (nM,uM,etc.)
                if allCytokinesHaveMWInDict:
                    cbaStandardsConcentrationMatrix[i,:] = np.multiply(cbaStandardsConcentrations[i,:],(concUnit/completeCytokineMWDict[cytokine]))
                    cbaStandardsConcentrationPlotPointsMatrix[i,:] = np.multiply(cbaStandardsConcentrationsPlotPoints[i,:],(concUnit/completeCytokineMWDict[cytokine]))
                else:
                    cbaStandardsConcentrationMatrix[i,:] = cbaStandardsConcentrations[i,:]
                    cbaStandardsConcentrationPlotPointsMatrix[i,:] = cbaStandardsConcentrationsPlotPoints[i,:]
                cbaStandardsMFIMatrix[i,:] = data[:,i]
                cbaStandardsMFIPlotPointsMatrix[i,:] = np.power(10,Hill(cbaStandardsConcentrationPlotPointsMatrix[i,:],*fittingParameters[i,:]))
                #Plot on log-log scale the experimental points and the curve fit line with previously determined curve fitting parameters
                #plt.loglog(cbaStandardsConcentrations,data[:,i],'o',color=color_list[i,:],label=listOfCytokines[i])
                #plt.loglog(cbaStandardsConcentrationsPlotPoints,np.power(10,Hill(convertedCBAStandardsPlotPoints,*fittingParameters[i,:])))
                #'_fit; R2 = '+str(rsquared)

                #Get LOD for each cytokine calibration curve (aka the linear range of the calibration curve)
                backgroundGFI = fittingParameters[i,3]
                amplitudeGFI = fittingParameters[i,0]
                
                #Approximate LOD by determining concentration values at LOD% and 1-LOD% (3% and 97%) of curve. Must be used on log10(curve), as calibration curve is plotted in logscale
                LODpercent = 0.03
                #LOD% more than background GFI used for lower LOD GFI
                lowerGFILOD = math.log10(10**((1+LODpercent)*math.log10(backgroundGFI)))
                #LOD% less than maximum GFI (Background + amplitude) used for upper LOD GFI
                upperGFILOD = math.log10(10**((1-LODpercent)*math.log10(amplitudeGFI+backgroundGFI)))
                #Log10(upper/lowerGFILOD) converted back into normal GFI by 10 to its power, then fed into inverse hill equation with current cytokine fitting parameters to obtain corresponding concentration values
                lowerConcLOD = InverseHill(lowerGFILOD,fittingParameters[i,:])
                upperConcLOD = InverseHill(upperGFILOD,fittingParameters[i,:])
                #Create dict with keys as cytokines, values as GFI/conc LODs
                concLOD[i,:] = np.array([10**lowerGFILOD,10**upperGFILOD,lowerConcLOD,upperConcLOD])
            flattenedMatrix = cbaStandardsMFIMatrix.flatten()
            reshapedMatrix = np.reshape(flattenedMatrix,(numberOfCalibrationSamples,len(cytokines)),order='F')
            flattenedMatrix2 = cbaStandardsMFIPlotPointsMatrix.flatten()
            reshapedMatrix2 = np.reshape(flattenedMatrix2,(numberOfPlotPoints,len(cytokines)),order='F')
            flattenedMatrix3 = cbaStandardsConcentrationMatrix.flatten()
            reshapedMatrix3 = np.reshape(flattenedMatrix3,(numberOfCalibrationSamples,len(cytokines)),order='F')
            flattenedMatrix4 = cbaStandardsConcentrationPlotPointsMatrix.flatten()
            reshapedMatrix4 = np.reshape(flattenedMatrix4,(numberOfPlotPoints,len(cytokines)),order='F')
            realCytokineList = cytokines
            dataValsList = []
            plotPointsList = []
            for j in range(1,numberOfCalibrationSamples+1):
                dataValsList.append([j])
            for j in range(1,numberOfPlotPoints+1):
                plotPointsList.append([j])
            dataValsIndex = pd.MultiIndex.from_tuples(dataValsList,names=['Standard'])
            plotPointsIndex = pd.MultiIndex.from_tuples(plotPointsList,names=['Standard'])
            currentCBAStandardsMFIDf = pd.DataFrame(reshapedMatrix,index=dataValsIndex,columns=realCytokineList)
            currentCBAPlotPointsMFIDf = pd.DataFrame(reshapedMatrix2,index=plotPointsIndex,columns=realCytokineList)
            currentCBAStandardsConcentrationDf = pd.DataFrame(reshapedMatrix3,index=dataValsIndex,columns=realCytokineList)
            currentCBAPlotPointsConcentrationDf = pd.DataFrame(reshapedMatrix4,index=plotPointsIndex,columns=realCytokineList)
           
            currentCBAStandardsMFIDf.columns.name = 'Cytokine' 
            currentCBAPlotPointsMFIDf.columns.name = 'Cytokine'
            currentCBAStandardsConcentrationDf.columns.name = 'Cytokine' 
            currentCBAPlotPointsConcentrationDf.columns.name = 'Cytokine'
            mic1 = pd.Index(cytokines, name='Cytokine')
            fittingParametersDf = pd.DataFrame(fittingParameters,index=mic1,columns=['Amplitude','EC50','HillCoeff','Background'])
            mic2 = pd.MultiIndex.from_tuples([['MFI','Lower'],['MFI','Upper'],['Concentration','Lower'],['Concentration','Upper']])
            concLODDf = pd.DataFrame(concLOD,index=mic1,columns=mic2)

            rsquaredDf = pd.DataFrame(rsquared_kit, index=mic1, columns=['R squared'])

            concLODList.append(concLODDf)
            fittingParametersList.append(fittingParametersDf)
            cbaStandardsMFIList.append(currentCBAStandardsMFIDf)
            cbaPlotPointsMFIList.append(currentCBAPlotPointsMFIDf)
            cbaStandardsConcentrationList.append(currentCBAStandardsConcentrationDf)
            cbaPlotPointsConcentrationList.append(currentCBAPlotPointsConcentrationDf)
            rsquaredList.append(rsquaredDf)

        #fullFittingParametersDf = pd.concat(fittingParametersList,keys=kitNames,names=['Kit Name'])
        fullConcLODDf = pd.concat(concLODList)
        fullFittingParametersDf = pd.concat(fittingParametersList)
        fullRsquaredDf = pd.concat(rsquaredList)
        #fullConcLODDf = pd.concat(concLODList)
        
        fullCBAStandardsMFIDf = pd.concat(cbaStandardsMFIList,keys=kitNames,names=['Kit Name'],axis=1)
        fullCBAPlotPointsMFIDf = pd.concat(cbaPlotPointsMFIList,keys=kitNames,names=['Kit Name'],axis=1)
        fullCBAStandardsConcentrationDf = pd.concat(cbaStandardsConcentrationList,keys=kitNames,names=['Kit Name'],axis=1)
        fullCBAPlotPointsConcentrationDf = pd.concat(cbaPlotPointsConcentrationList,keys=kitNames,names=['Kit Name'],axis=1)
        fullCBAStandardsList = [fullCBAStandardsMFIDf.stack().stack(),fullCBAStandardsConcentrationDf.stack().stack()]
        fullCBAPlotPointsList = [fullCBAPlotPointsMFIDf.stack().stack(),fullCBAPlotPointsConcentrationDf.stack().stack()]
        fullCBAStandardsDf = pd.concat(fullCBAStandardsList,axis=1,keys=[yaxistitle,xaxistitle])
        fullCBAPlotPointsDf = pd.concat(fullCBAPlotPointsList,axis=1,keys=[yaxistitle,xaxistitle])
         
        plottingPointsDf = fullCBAPlotPointsDf.reset_index()
        plottingStandardsDf = fullCBAStandardsDf.reset_index()
        
        numCyt = len(pd.unique(plottingPointsDf['Cytokine']))
        if numCyt <= 13:
            fullpalette = sns.color_palette(sns.color_palette(cc.glasbey),numCyt)
            g = sns.relplot(data=plottingPointsDf,x=xaxistitle,y=yaxistitle,hue='Cytokine',col='Kit Name',kind='line',col_order=pd.unique(plottingPointsDf['Kit Name']),hue_order=pd.unique(plottingPointsDf['Cytokine']),height=7,palette=fullpalette)
            #Plot vertical lines at lower and upper concentration limits of detection
            colorDict = {}
            for j,cytokine in enumerate(pd.unique(plottingPointsDf['Cytokine'])):
                colorDict[cytokine] = fullpalette[j]
            for axis,kitName in zip(g.axes.flat,pd.unique(plottingPointsDf['Kit Name'])):
                currentpalette = []
                for cytokine in pd.unique(plottingStandardsDf.query("`Kit Name` == @kitName")['Cytokine']):
                    currentColor = colorDict[cytokine]
                    currentpalette.append(currentColor)
                    cytokineLODValues = fullConcLODDf.loc[cytokine,:]['Concentration']
                    axis.axvline(x=cytokineLODValues['Lower'],color=currentColor,linestyle=':')
                    axis.axvline(x=cytokineLODValues['Upper'],color=currentColor,linestyle=':')
                g2 = sns.scatterplot(data=plottingStandardsDf[plottingStandardsDf['Kit Name'] == kitName],x=xaxistitle,y=yaxistitle,hue='Cytokine',ax=axis,legend=False,palette=currentpalette)
                axis.set_xscale('log')
                axis.set_yscale('log')
                if not allCytokinesHaveMWInDict:
                    axis.set_ylabel('Concentration of Cytokine Standards (pg/mL)')
                    concUnitPrefix = 'pg_mL'
            plt.savefig('plots/calibrationCurves-'+folderName+'-'+concUnitPrefix+'.png',bbox_inches='tight', dpi=200)
        else:
            fullpalette = sns.color_palette(sns.color_palette(),len(pd.unique(plottingPointsDf['Kit Name'])))
            colorDict = {}
            g = sns.relplot(data=plottingPointsDf,x=xaxistitle,y=yaxistitle,hue='Kit Name',col='Cytokine',kind='line',hue_order=pd.unique(plottingPointsDf['Kit Name']),col_order=pd.unique(plottingPointsDf['Cytokine']),col_wrap=math.floor(np.sqrt(numCyt)),height=12,palette=fullpalette,facet_kws={'sharey':False, 'sharex':False})
            for j,kit in enumerate(pd.unique(plottingPointsDf['Kit Name'])):
                    colorDict[kit] = fullpalette[j]
            for axis,cytName in zip(g.axes.flat,pd.unique(plottingPointsDf['Cytokine'])):
                currentpalette = []
                for kit in pd.unique(plottingStandardsDf.query("Cytokine == @cytName")['Kit Name']):
                    currentColor = colorDict[kit]
                    currentpalette.append(currentColor)
                    cytokineLODValues = fullConcLODDf.loc[cytName,:]['Concentration']
                    axis.axvline(x=cytokineLODValues['Lower'],color=currentColor,linestyle=':')
                    axis.axvline(x=cytokineLODValues['Upper'],color=currentColor,linestyle=':')
                g2 = sns.scatterplot(data=plottingStandardsDf[plottingStandardsDf['Cytokine'] == cytName],x=xaxistitle,y=yaxistitle,hue='Kit Name',ax=axis,legend=False)
                axis.set_xscale('log')
                axis.set_yscale('log')
                if not allCytokinesHaveMWInDict:
                    axis.set_ylabel('Concentration of Cytokine Standards (pg/mL)')
                    concUnitPrefix = 'pg_mL'
            plt.tight_layout()
            g.fig.savefig('plots'+dirSep+'calibrationCurves-'+folderName+'-'+concUnitPrefix+'.png',bbox_inches='tight', dpi=200)
        #Save fitting parameters and LOD for curve fit for each cytokine
        with open('misc'+dirSep+'fittingParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "wb") as f:
            pickle.dump(fullFittingParametersDf, f)
        with open('misc'+dirSep+'LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "wb") as f:
            pickle.dump(fullConcLODDf, f)
        with open('misc/rsquaredValues-'+folderName+'-'+concUnitPrefix+'.pkl', "wb") as f:
            pickle.dump(fullRsquaredDf, f)
    return unparsedCytokines

def createCytokineDataFrame(folderName,finalDataFrame,concUnitPrefix):
         
        columnName = finalDataFrame.columns.name
        with open('outputData'+dirSep+'pickleFiles'+dirSep+'cytokineGFIPickleFile-'+folderName+'.pkl', "wb") as f:
            pickle.dump(finalDataFrame, f)
        
        fittingParameters = pickle.load(open('misc'+dirSep+'fittingParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
        LODParameters = pickle.load(open('misc'+dirSep+'LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
        #Begin converting GFI dataframe into corresponding concentration dataframe
        concentrationList = []
        #Step through dataframe one cytokine at a time
        idx = pd.IndexSlice
        for cytokine in pd.unique(finalDataFrame.index.get_level_values(0)):
            #Retrieve LODs for current cytokine (from constructed calibration curve)
            lowerGFILOD = LODParameters.loc[cytokine,idx['MFI','Lower']]
            upperGFILOD = LODParameters.loc[cytokine,idx['MFI','Upper']]
            lowerConcLOD = LODParameters.loc[cytokine,idx['Concentration','Lower']]
            upperConcLOD = LODParameters.loc[cytokine,idx['Concentration','Upper']]
            smallConcentrationMatrix = np.zeros(finalDataFrame.loc[cytokine].shape)
            #Loop through every value in current cytokine's portion of the dataframe
            for i in range(0,finalDataFrame.loc[cytokine].values.shape[0]):
                for j in range(0,finalDataFrame.loc[cytokine].values.shape[1]):
                    currentGFIval = finalDataFrame.loc[cytokine].values[i,j]
                    if currentGFIval > upperGFILOD: #If intensity is greater than upper GFI LOD
                        currentConcVal = upperConcLOD #Concentration is equal to upper concentration LOD
                    elif currentGFIval <= upperGFILOD and currentGFIval >= lowerGFILOD: #if intensity is between upper and lower GFI LODs
                        currentConcVal = InverseHill(np.log10(currentGFIval),fittingParameters.loc[cytokine,:].values) #Use previous hill fit parameters for the cytokine to obtain concentration
                    else: #If intensity is less than background GFI LOD
                        currentConcVal = lowerConcLOD #Concentration is equal to lower concentration LOD
                    smallConcentrationMatrix[i,j] = currentConcVal
            concentrationList.append(smallConcentrationMatrix)

        concentrationMatrix = np.vstack(concentrationList)
        finalDataFrameConcentration = pd.DataFrame(concentrationMatrix,index=finalDataFrame.index,columns=finalDataFrame.columns)
        finalDataFrameConcentration.columns.name = columnName 
        return finalDataFrameConcentration
