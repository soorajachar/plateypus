#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,glob
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import numpy as np
import pandas as pd
import seaborn as sns
import tkinter as tk
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from itertools import groupby
from plateypus.dataprocessing.miscFunctions import Hill,InverseHill,r_squared,cleanUpFlowjoCSV
idx = pd.IndexSlice

#Th1/2/17 Mouse  BD Biosciences CBA Kit Cytokines
bdMouseThKitDict = {'IFNg':17200,'IL-2':17200,'IL-4':14000,'IL-6':21900,'IL-10':18900,'IL-17A':15500,'TNFa':17500}
#Inflammatory Mouse  BD Biosciences CBA Kit Cytokines
bdMouseInfKitDict = {'IFNg':17200,'IL-6':21900,'IL-10':18900,'TNFa':17500,'IL-12p70':57480,'MCP-1':16000}
#Inflammatory Human BD Biosciences CBA Kit Cytokines
bdHumanInfKitDict = {'IL-12p70':70000,'IL-1B':17300,'IL-8':8904}
#Chemokine Human BD Biosciences CBA Kit Cytokines
bdHumanChemokineKitDict = {'IP-10':86000,'MCP-1':13000,'MIG':11700,'RANTES':7809,'IL-8':8904}
#Legendplex 13-plex mouse th kit
legendPlexMouseThKitDict = {'IFNg':15652,'IL-5':26200,'TNFa':25896,'IL-2':17231,'IL-6':21709,'IL-4':13500,'IL-10':20641,'IL-9':14300,'IL-17A':14978,'IL-17F':14900,'IL-21':14400,'IL-22':16800,'IL-13':12300}
#Legendplex 13-plex mouse macrophage kit
legendPlexMouseMacrophageKitDict = {'G-CSF':18800,'GM-CSF':14000,'IL-12':75000,'IL-12p40':40000,'IL-18':18000,'IL-1b':17500,'IL-23':55000,'IL-34':39000,'IL-7':17000,'KC':11000,'MCP-1':12000,'TARC':8000,'TGFb':25000}

#Combine all cytokine bead kit molar masses:
listOfKitDicts = [legendPlexMouseThKitDict,bdMouseThKitDict,bdHumanInfKitDict,bdMouseInfKitDict,bdHumanChemokineKitDict,legendPlexMouseMacrophageKitDict]
completeCytokineMWDict = {}
for dictName in listOfKitDicts:
    completeCytokineMWDict.update(dictName)
"""
#Standard BD Biosciences Th1/2/17 Kit Cytokines
listOfCytokines1=['IFNg','IL-2','IL-4','IL-6','IL-10','IL-17A','TNFa']
MWofCytokines1=[17200,17200,14000,21900,18900,15500,17500] #g/mol, correct one

#Standard BD Biosciences Inflamamatory Kit Cytokines
listOfCytokines2=['IL-12p70','IL-1B','IL-8']
MWofCytokines2=[70000,17300,8904] #g/mol, correct one

#Standard Mouse Th Legendplex CBA Kit Cytokines
#Up to IL-4 is A (lower FSC and SSC), after that is B
legendPlexListOfCytokines1=['IFNg','IL-5','TNFa','IL-2','IL-6','IL-4','IL-10','IL-9','IL-17A','IL-17F','IL-21','IL-22','IL-13']
"""

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
         
        def collectInputs():
            numCalibrationSamples = int(t1.get())
            initialStandardVolume = float(t2.get())
            calibrationParameterDict = {'Number': numCalibrationSamples,'Volume': initialStandardVolume}
            with open('misc/CBAcalibrationParameters-'+folderName+'.json','w') as f:
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
    with open('inputData/bulkCSVFiles/'+fileName, 'r') as istr:
        with open('inputData/bulkCSVFiles/'+fileName, 'w') as ostr:
            for line in istr:
                if line[-1] != ',':
                    line = line.rstrip('\n') + ','
                    print(line, file=ostr)

#Create Calibration Curves, obtain LODs (limits of detection) of experiment
def calibrateExperiment(folderName,secondPath,concUnit,concUnitPrefix,numberOfCalibrationSamples,initialStandardVolume):
    #Get cytokine calibration curve data
    tempExperimentParameters = {'overallPlateDimensions':[8,12]}
    calibrationFileNames = glob.glob('inputData/bulkCSVFiles/Calibration*')
    calibrationNames = []
    kitNames = []
    for calibrationFileName in calibrationFileNames:
        #performCommaCheck(calibrationFileName.split('/')[-1])
        newName = calibrationFileName.split('.')[0].split('_')[0].split('/')[-1]
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
        
    allCytokinesHaveMWInDict = True 
    for calibration in sortedData:
        cytokines = parseCytokineCSVHeaders(calibration.columns)
        for cytokine in cytokines:
            if cytokine[0] not in completeCytokineMWDict:
                allCytokinesHaveMWInDict = False

    for calibration in sortedData:
        data = np.array(calibration.values[:,1:-1],dtype=float)
        cytokines = parseCytokineCSVHeaders(calibration.columns)
        fittingParameters = np.zeros((data.shape[1],4))
        concLOD = np.zeros((data.shape[1],4))
        #Initial concentration all cytokine standards is given by CBA kit manual as 5000 pGg/mL: when standards are diluted in 2mL
        conc = 5000 #pg/mL
        serialDilutionFactor = 2 #1:serialDilutionFactor dilution between each standard well
        #Smaller initial dilution (0.5mL instead of 2mL for example) increase the initial concentration of the first calibration sample
        initialConc = (conc*1e-12) /((initialStandardVolume*1e-3)/2) #g/L (pg/mL * 1e-12 g/pg)/(1e-3 L/mL)
        #Calibration samples are always diluted by a factor of serialdilutionFactor (so with 12 calibration samples, the last sample is (serialDilutionFactor^-11) the concentration of the first, which is pure standard (2^0)
        cbaStandardsConcentrations = np.flipud(initialConc*np.power(serialDilutionFactor,np.linspace(-numberOfCalibrationSamples+1,0,numberOfCalibrationSamples)))
        #More x values along the above concentration bounds are sampled to use to construct calibration curve. Plot points are extended slightly at high range to allow visualization of upper LOD not accessible with experimental dilution
        cbaStandardsConcentrationsPlotPoints = np.flipud(initialConc*np.power(2,np.linspace(-numberOfCalibrationSamples+1,4,numberOfPlotPoints)))
        
        cbaStandardsConcentrationMatrix = np.zeros([len(cytokines),cbaStandardsConcentrations.shape[0]]) 
        cbaStandardsConcentrationPlotPointsMatrix = np.zeros([len(cytokines),cbaStandardsConcentrationsPlotPoints.shape[0]])
        cbaStandardsMFIMatrix = np.zeros([len(cytokines),cbaStandardsConcentrations.shape[0]]) 
        cbaStandardsMFIPlotPointsMatrix = np.zeros([len(cytokines),cbaStandardsConcentrationsPlotPoints.shape[0]])
        color_list = sns.color_palette(sns.color_palette(),len(cytokines))
        for i,cytokineList in enumerate(cytokines):
            cytokine = cytokineList[0]
            #amplitude bounded from range/2 to range*2, EC50 bounded from minimum to maximum standard concentration tested, Hill coefficient bounded from 0 to 2, Background bounded from 0 to minimum GFI*2
            lowerCurveFitBounds = [(np.max(data[:,i])-np.min(data[:,i]))/2,np.min(cbaStandardsConcentrations),0,0]
            upperCurveFitBounds = [(np.max(data[:,i])-np.min(data[:,i]))*2, np.max(cbaStandardsConcentrations), 2,np.min(data[:,i])*2]
            #use scipy curve fit to determine best hill equation fit for data, searching within the bounds given above
            popt,pcov = curve_fit(Hill, cbaStandardsConcentrations,np.log10(data[:,i]),sigma=np.log10(data[:,i]),bounds=(lowerCurveFitBounds,upperCurveFitBounds))
            rsquared = round(r_squared(cbaStandardsConcentrations,np.log10(data[:,i]),Hill,popt),3)
            rsquaredList.append(rsquared)
            for j in range(len(popt)):  
                #Convert just ec50 value to desired units (nM,uM etc) if cytokine has a molar mass in dict
                if j == 1 and allCytokinesHaveMWInDict:
                    fittingParameters[i,j] = np.multiply(popt[j],(concUnit/completeCytokineMWDict[cytokine]))
                #other values in 4 parameter logistic equation are tied to intensity y-value, which doesn't change, or are the hill coefficient, which is completely separate, so parameters are kept the same
                else:
                    fittingParameters[i,j] = popt[j]
            
            #Convert x values of experimental data points and curve fit points to desired units (nM,uM,etc.)
            if allCytokinesHaveMWInDict:
                cbaStandardsConcentrationMatrix[i,:] = np.multiply(cbaStandardsConcentrations,(concUnit/completeCytokineMWDict[cytokine]))
                cbaStandardsConcentrationPlotPointsMatrix[i,:] = np.multiply(cbaStandardsConcentrationsPlotPoints,(concUnit/completeCytokineMWDict[cytokine]))
            else:
                cbaStandardsConcentrationMatrix[i,:] = cbaStandardsConcentrations
                cbaStandardsConcentrationPlotPointsMatrix[i,:] = cbaStandardsConcentrationsPlotPoints
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
        realCytokineList = []
        for cytokine in cytokines:
            realCytokineList.append(cytokine[0])
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
        mic1 = pd.MultiIndex.from_tuples(cytokines,names=['Cytokine'])
        fittingParametersDf = pd.DataFrame(fittingParameters,index=mic1,columns=['Amplitude','EC50','HillCoeff','Background'])
        mic2 = pd.MultiIndex.from_tuples([['MFI','Lower'],['MFI','Upper'],['Concentration','Lower'],['Concentration','Upper']])
        concLODDf = pd.DataFrame(concLOD,index=mic1,columns=mic2)

        concLODList.append(concLODDf)
        fittingParametersList.append(fittingParametersDf)
        cbaStandardsMFIList.append(currentCBAStandardsMFIDf)
        cbaPlotPointsMFIList.append(currentCBAPlotPointsMFIDf)
        cbaStandardsConcentrationList.append(currentCBAStandardsConcentrationDf)
        cbaPlotPointsConcentrationList.append(currentCBAPlotPointsConcentrationDf)

    #fullFittingParametersDf = pd.concat(fittingParametersList,keys=kitNames,names=['Kit Name'])
    #fullConcLODDf = pd.concat(concLODList,keys=kitNames,names=['Kit Name'])
    fullFittingParametersDf = pd.concat(fittingParametersList)
    fullConcLODDf = pd.concat(concLODList)
    
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
    if numCyt <= 10:
        fullpalette = sns.color_palette(sns.color_palette(),numCyt)
    else:
        fullpalette = sns.color_palette('hls',numCyt)
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
            axis.axvline(x=cytokineLODValues['Lower'].values,color=currentColor,linestyle=':')
            axis.axvline(x=cytokineLODValues['Upper'].values,color=currentColor,linestyle=':')
        g2 = sns.scatterplot(data=plottingStandardsDf[plottingStandardsDf['Kit Name'] == kitName],x=xaxistitle,y=yaxistitle,hue='Cytokine',ax=axis,legend=False,palette=currentpalette)
        axis.set_xscale('log')
        axis.set_yscale('log')
    plt.savefig('plots/calibrationCurves-'+folderName+'-'+concUnitPrefix+'.png',bbox_inches='tight')
    #Save fitting parameters and LOD for curve fit for each cytokine
    with open('misc/fittingParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "wb") as f:
        pickle.dump(fullFittingParametersDf, f)
    with open('misc/LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "wb") as f:
        pickle.dump(fullConcLODDf, f)

def createCytokineDataFrame(folderName,finalDataFrame,concUnitPrefix):
         
        columnName = finalDataFrame.columns.name
        with open('outputData/pickleFiles/cytokineGFIPickleFile-'+folderName+'.pkl', "wb") as f:
            pickle.dump(finalDataFrame, f)
        
        fittingParameters = pickle.load(open('misc/fittingParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
        LODParameters = pickle.load(open('misc/LODParameters-'+folderName+'-'+concUnitPrefix+'.pkl', "rb"))
        #Begin converting GFI dataframe into corresponding concentration dataframe
        concentrationList = []
        #Step through dataframe one cytokine at a time
        for cytokine in pd.unique(finalDataFrame.index.get_level_values(0)):
            #Retrieve LODs for current cytokine (from constructed calibration curve)
            lowerGFILOD = LODParameters.loc[cytokine,idx['MFI','Lower']].values[0]
            upperGFILOD = LODParameters.loc[cytokine,idx['MFI','Upper']].values[0] 
            lowerConcLOD = LODParameters.loc[cytokine,idx['Concentration','Lower']].values[0]
            upperConcLOD = LODParameters.loc[cytokine,idx['Concentration','Upper']].values[0] 
            smallConcentrationMatrix = np.zeros(finalDataFrame.loc[cytokine].shape)
            #Loop through every value in current cytokine's portion of the dataframe
            for i in range(0,finalDataFrame.loc[cytokine].values.shape[0]):
                for j in range(0,finalDataFrame.loc[cytokine].values.shape[1]):
                    currentGFIval = finalDataFrame.loc[cytokine].values[i,j]
                    if currentGFIval > upperGFILOD: #If intensity is greater than upper GFI LOD
                        currentConcVal = upperConcLOD #Concentration is equal to upper concentration LOD
                    elif currentGFIval <= upperGFILOD and currentGFIval >= lowerGFILOD: #if intensity is between upper and lower GFI LODs
                        currentConcVal = InverseHill(np.log10(currentGFIval),*fittingParameters.loc[cytokine,:].values) #Use previous hill fit parameters for the cytokine to obtain concentration
                    else: #If intensity is less than background GFI LOD
                        currentConcVal = lowerConcLOD #Concentration is equal to lower concentration LOD
                    smallConcentrationMatrix[i,j] = currentConcVal
            concentrationList.append(smallConcentrationMatrix)

        concentrationMatrix = np.vstack(concentrationList)
        finalDataFrameConcentration = pd.DataFrame(concentrationMatrix,index=finalDataFrame.index,columns=finalDataFrame.columns)
        finalDataFrameConcentration.columns.name = columnName 
        return finalDataFrameConcentration
