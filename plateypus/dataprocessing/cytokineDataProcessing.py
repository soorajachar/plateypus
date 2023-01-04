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
import os
if os.name == 'nt':
    dirSep = '\\'
else:
    dirSep = '/'
from plateypus.dataprocessing.miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV, setMaxWidth
idx = pd.IndexSlice

#Th1/2/17 Mouse  BD Biosciences CBA Kit Cytokines
bdMouseThKitDict = {'IFNg':17200,'IL-2':17200,'IL-4':14000,'IL-6':21900,'IL-10':18900,'IL-17A':15500,'TNFa':17500}
#Inflammatory Mouse  BD Biosciences CBA Kit Cytokines
bdMouseInfKitDict = {'IFNg':17200,'IL-6':21900,'IL-10':18900,'TNFa':17500,'IL-12p70':57480,'MCP-1':16000}
#Th1/2/17 Human  BD Biosciences CBA Kit Cytokines
bdHumanThKitDict = {'IFNg':16900,'IL-2':15386,'IL-4':15000,'IL-6':21000,'IL-10':18600,'IL-17A':30700,'TNFa':17500}
#Inflammatory Human BD Biosciences CBA Kit Cytokines
bdHumanInfKitDict = {'IL-12p70':70000,'IL-1B':17300,'IL-8':8904}
#Chemokine Human BD Biosciences CBA Kit Cytokines
bdHumanChemokineKitDict = {'IP-10':8600,'MCP-1':13000,'MIG':11700,'RANTES':7809,'IL-8':8904}
#Legendplex 13-plex mouse th kit
legendPlexMouseThKitDict = {'IFNg':15652,'IL-5':26200,'TNFa':25896,'IL-2':17231,'IL-6':21709,'IL-4':13500,'IL-10':20641,'IL-9':14300,'IL-17A':14978,'IL-17F':14900,'IL-21':14400,'IL-22':16800,'IL-13':12300}
#Legendplex 13-plex mouse macrophage kit
legendPlexMouseMacrophageKitDict = {'G-CSF':18800,'GM-CSF':14000,'IL-12':75000,'IL-12p40':40000,'IL-18':18000,'IL-1b':17500,'IL-23':55000,'IL-34':39000,'IL-7':17000,'KC':11000,'MCP-1':12000,'TARC':8000,'TGFb':25000}
# Human Soluble Protein Flex CBA Kit
bdHumanFlexKit = {'Angiogenin': 14000,'CD121a': 37400,'CD121b': 38800,'CD178': 18000,'CD40L': 16300,'CD54': 95000,'CD62L': 60190,'Eotaxin': 8400,'FGF': 16310,'Fractalkine': 8770,'G-CSF': 18600,\
    'GM-CSF': 14500,'Granzyme A': 28000,'Granzyme B': 27500,'IFNg': 16900,'IL-10': 18600,'IL-11': 19100,'IL-12': 57000,'IL-17A': 30700,'IL-17F': 14903,'IL-1A': 18047,'IL-1B': 17300,'IL-2': 15386,\
    'IL-21': 15500,'IL-3': 15000,'IL-4': 15000,'IL-5': 26522,'IL-6': 21000,'IL-7': 17400,'IL-8': 8904,'IL-9': 14000,'IP-10': 8600,'LT-Alp': 18600,'MIG': 11700,'MIP-1A': 7500,'RANTES': 7809,'TGF-B': 12700,'TNF': 17500,'TNFRI': 22700,'TNFRII': 26600,'VEGF': 19000}
# Mouse Soluble Protein Flex CBA Kit
bdMouseFlexKit = {'IFNg':17200, 'IL-2':17200, 'IL-5':26200, 'IL-4':1400, 'IL-3':15200, 'KC':11000, 'IL-6':21900, 'IL-21':14400, 'MCP-1':16000, 'IL-13':12300, 'IL-10':18900, 'IL-17A':15500, 'MIP-1a':7820, \
    'TNF':17500, 'IL-17F':14900, 'IL-12_IL-23p40':70000, 'RANTES':7800, 'MIG':12200, 'IL-1a':17500, 'IL-1b':17500, 'CD62E':115000, 'CD62L':76000}

listOfHumanKitDicts = [bdHumanThKitDict, bdHumanInfKitDict,bdHumanChemokineKitDict,bdHumanFlexKit]
listofMouseKitDicts = [legendPlexMouseThKitDict,bdMouseThKitDict,bdMouseInfKitDict,legendPlexMouseMacrophageKitDict, bdMouseFlexKit]

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

        l2 = tk.Label(mainWindow,text='Volume of initial CBA calibration solution: ').grid(row=0,column=1)
        t2 = tk.Entry(mainWindow)
        t2.grid(row=1,column=1,sticky=tk.W)

        l3 = tk.Label(mainWindow, text='Species of CBA samples: ').grid(row=0, column=3)
        species = ['Mouse', 'Human']
        speciesVar = tk.StringVar()
        t3 = tk.OptionMenu(mainWindow,speciesVar,*species)
        setMaxWidth(species, t3)
        t3.grid(row=1, column=3)
         
        def collectInputs():
            numCalibrationSamples = int(t1.get())
            initialStandardVolume = float(t2.get())
            speciesSamples = speciesVar.get()
            calibrationParameterDict = {'Number': numCalibrationSamples,'Volume': initialStandardVolume, 'Species':speciesSamples}
            with open('misc'+dirSep+'CBAcalibrationParameters-'+folderName+'.json','w') as f:
                json.dump(calibrationParameterDict,f)
            master.switch_frame(secondaryhomepage,folderName,expNum,ex_data,bPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

def parseCytokineCSVHeaders(columns):
    #,Beads/IFNg | Geometric Mean (YG586-A),Beads/IL-2 | Geometric Mean (YG586-A),Beads/IL-4 | Geometric Mean (YG586-A),Beads/IL-6 | Geometric Mean (YG586-A),Beads/IL-10 | Geometric Mean (YG586-A),Beads/IL-17A | Geometric Mean (YG586-A),Beads/TNFa | Geometric Mean (YG586-A),
    newMultiIndexList = []
    for column in columns[1:-1]:
        populationNameVsStatisticSplit = column.split(' | ')
        if 'Barcode' not in populationNameVsStatisticSplit[0]:
            cytokine = populationNameVsStatisticSplit[0].split('/')[-1]
        else:
            cytokine = populationNameVsStatisticSplit[0].split('/')[-2]
        newMultiIndexList.append([cytokine])
    return newMultiIndexList

def performCommaCheck(fileName):
    with open('inputData'+dirSep+'bulkCSVFiles'+dirSep+fileName, 'r') as istr:
        with open('inputData'+dirSep+'bulkCSVFiles'+dirSep+fileName, 'w') as ostr:
            for line in istr:
                if line[-1] != ',':
                    line = line.rstrip('\n') + ','
                    print(line, file=ostr)

#Create Calibration Curves, obtain LODs (limits of detection) of experiment
def calibrateExperiment(folderName,secondPath,concUnit,concUnitPrefix,numberOfCalibrationSamples,initialStandardVolume, species):
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
    completeCytokineMWDict = {}
    if species == 'Human':
        for dictName in listOfHumanKitDicts:
            completeCytokineMWDict.update(dictName)
    else:
        for dictName in listofMouseKitDicts:
            completeCytokineMWDict.update(dictName)
        
    allCytokinesHaveMWInDict = True 
    for calibration in sortedData:
        cytokines = parseCytokineCSVHeaders(calibration.columns)
        cytokines = [x[0] for x in cytokines]
        for cytokine in cytokines:
            if cytokine not in completeCytokineMWDict:
                allCytokinesHaveMWInDict = False

    for calibration in sortedData:
        data = np.array(calibration.values[:,1:-1],dtype=float)
        cytokines = parseCytokineCSVHeaders(calibration.columns)
        cytokines = [x[0] for x in cytokines]
        fittingParameters = np.zeros((data.shape[1],4))
        concLOD = np.zeros((data.shape[1],4))
        serialDilutionFactor = 2 # serialDilutionFactor dilution between each standard well

        
        if len(cytokines) > 12:
            if species == 'Human':
            #Initial concentration of cytokine standards given by individual Flex kit manuals in pg/mL when diluted in 4mL
                all_conc = {'Angiogenin': 2500,'CD121a': 10000,'CD121b': 10000,'CD178': 2500,'CD40L': 2500,'CD54': 10000,'CD62L': 10000,'Eotaxin': 2500,'FGF': 2500,'Fractalkine': 10000,'G-CSF': 2500,'GM-CSF': 2500,\
                    'Granzyme A': 10000,'Granzyme B': 10000,'IFNg': 2500,'IL-10': 2500,'IL-11': 10000,'IL-12': 10000,'IL-17A': 2500,'IL-17F': 2500,'IL-1A': 2500,'IL-1B': 2500,'IL-2': 2500,'IL-21': 10000,'IL-3': 2500,\
                    'IL-4': 2500,'IL-5': 2500,'IL-6': 2500,'IL-7': 2500,'IL-8': 2500,'IL-9': 2500,'IP-10': 2500,'LT-Alp': 2500,'MIG': 2500,'MIP-1A': 2500,'RANTES': 2500,'TGF-B': 10000,'TNF': 2500,'TNFRI': 10000,'TNFRII': 2500,'VEGF': 2500}
                conc = np.array([all_conc[cyt] for cyt in cytokines])
                #Smaller initial dilution (2mL instead of 4mL for example) increase the initial concentration of the first calibration sample
                initialConc = (conc*1e-12) /((initialStandardVolume*1e-3)/4)
            else:
                all_conc = {'IFNg':2500, 'IL-2':2500, 'IL-5':2500, 'IL-4':2500, 'IL-3':2500, 'KC':2500, 'IL-6':2500, 'IL-21':10000, 'MCP-1':2500, 'IL-13':2500, 'IL-10':2500, 'IL-17A':2500, 'MIP-1a':2500, \
                    'TNF':2500, 'IL-17F':2500, 'IL-12_IL-23p40':10000, 'RANTES':2500, 'MIG':2500, 'IL-1a':2500, 'IL-1b':2500, 'CD62E':10000, 'CD62L':10000}
                conc = np.array([all_conc[cyt] for cyt in cytokines])
                #Smaller initial dilution (2mL instead of 4mL for example) increase the initial concentration of the first calibration sample
                initialConc = (conc*1e-12) /((initialStandardVolume*1e-3)/4)
        else:
            #Initial concentration all cytokine standards is given by CBA kit manual as 5000 pg/mL: when standards are diluted in 2mL
            conc = np.tile([5000], len(cytokines))
            #Smaller initial dilution (0.5mL instead of 2mL for example) increase the initial concentration of the first calibration sample
            initialConc = (conc*1e-12) /((initialStandardVolume*1e-3)/2) #g/L (pg/mL * 1e-12 g/pg)/(1e-3 L/mL)
        
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
    if numCyt <= 12:
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
