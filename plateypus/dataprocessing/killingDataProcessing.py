#!/usr/bin/env python3 
import string,os
import pandas as pd
import numpy as np
import os
if os.name == 'nt':
    dirSep = '\\'
else:
    dirSep = '/'
from plateypus.dataprocessing.miscFunctions import reorderDfByInputOrder

def generateBulkKillingStatistics(experimentParameters,folderName,experimentNumber,dataType,layoutDict):
        numRowPlates = 1 
        numColumnPlates = experimentParameters['numPlates']

        #Combine plate and well IDs into a single ID val for every single sample in the experiment
        identificationMatrix = np.empty(layoutDict['plateID'].shape,dtype=object)
        for row in range(identificationMatrix.shape[0]):
            for col in range(identificationMatrix.shape[1]):
                wellID = layoutDict['wellID'][row,col]
                plateID = layoutDict['plateID'][row,col]
                fullID = plateID+'-'+wellID
                identificationMatrix[row,col] = fullID

        plateNames = np.unique(layoutDict['plateID'])
        
        levelLabelDict = experimentParameters['levelLabelDict']
        if 'Time' in levelLabelDict:
            del levelLabelDict['Time']

        plateLayoutDict = layoutDict['keys']
        blankLayout = layoutDict['blank'].astype(object)
        
        fullExperimentDf,plateLayoutDf = createIncucyteDataframe(plateNames,{**plateLayoutDict,**{'Blank':blankLayout}},{**levelLabelDict,**{'Blank':['Yes','No']}}) 
        
        fullExperimentDf = reorderDfByInputOrder(experimentParameters,fullExperimentDf)
        return fullExperimentDf

def createIncucyteDataframe(plateNames,plateLayoutDict,levelLabelDict,unwrappingOrder='col'):
    path = 'inputData'+dirSep+'bulkCSVFiles'+dirSep
    if 'A1_killing.csv' in os.listdir(path):
        csvBool = True
    else:
        csvBool = False
    #Create incucyte data matrix
    dfList = []
    for plateName in plateNames:
        if csvBool:
            plateDf = pd.read_csv(path+plateName+'_killing.csv',header=1).iloc[:,2:]
        else:
            plateDf = pd.read_excel(path+plateName+'_killing.xlsx',header=1).iloc[:,2:]
        dfList.append(plateDf.dropna())
    
    inputDf = pd.concat(dfList,axis=1,keys=string.ascii_uppercase[:len(dfList)]).dropna()
    if csvBool:
        inputDf.index = pd.read_csv(path+plateNames[0]+'_killing.csv',header=1).iloc[:,1].dropna().values
    else:
        inputDf.index = pd.read_excel(path+plateNames[0]+'_killing.xlsx',header=1).iloc[:,1].dropna().values
    
    inputDf.index.name = 'Time'
    
    unwrappingOrderDict = {'col':'F','row':'C'}
    flatteningOrder = unwrappingOrderDict[unwrappingOrder]

    inputPlateLayout = {}
    for levelKey,level in zip(plateLayoutDict,levelLabelDict):
      inputPlateLayout[level] = plateLayoutDict[levelKey]
    
    # Reshape lists so that they are in row/column major format
    reshapedLevelValues = []
    for level in inputPlateLayout:
        plateLayout = inputPlateLayout[level]
        levelValues = levelLabelDict[level]
        for row in range(plateLayout.shape[0]):
            for col in range(plateLayout.shape[1]):
                if level != 'Blank':
                    if plateLayout[row,col] == 'blank':
                        plateLayout[row,col] = 'blank'
                    else:
                        plateLayout[row,col] = levelValues[plateLayout[row,col]]
                else:
                    if plateLayout[row,col] == 0:
                        plateLayout[row,col] = 'blank'
                    else:
                        plateLayout[row,col] = levelValues[plateLayout[row,col]]
        newLayout = [x for x in plateLayout.flatten(order=flatteningOrder) if x != 'blank']
        reshapedLevelValues.append(newLayout)

    #Create multiIndex
    multiIndex = pd.MultiIndex.from_arrays(reshapedLevelValues,names=inputPlateLayout.keys())
    #Drop unused wells (in block plates) from multiIndex
    trueMultiIndex = pd.MultiIndex.from_tuples([x for x in multiIndex.values.tolist() if '' not in x],names=inputPlateLayout.keys())

    #Create plate layout dataframe
    plateLayoutDf = pd.DataFrame(trueMultiIndex.values.tolist(),index=inputDf.columns,columns=list(inputPlateLayout.keys()))
    plateLayoutDf.index.names = ['Plate','Well']

    #Create final dataframe, drop blanks, and swap columns/rows
    if 'Blank' in list(inputPlateLayout.keys()):
        if 'Yes' in inputPlateLayout['Blank'].flatten():
            incucyteDf = pd.DataFrame(inputDf.values,index=inputDf.index,columns=trueMultiIndex).drop('Yes',level='Blank',axis=1).droplevel('Blank',axis=1).replace('na',np.nan)#.astype(float).T
        else:
            incucyteDf = pd.DataFrame(inputDf.values,index=inputDf.index,columns=trueMultiIndex).droplevel('Blank',axis=1).replace('na',np.nan)#.astype(float).T
    else:
        incucyteDf = pd.DataFrame(inputDf.values,index=inputDf.index,columns=trueMultiIndex).replace('na',np.nan)#.astype(float).T
    
    #Correct missing values to zero
    for i in range(incucyteDf.shape[0]):
        for j in range(incucyteDf.shape[1]):
            if type(incucyteDf.iloc[i,j]) == str:
                incucyteDf.iloc[i,j] = 0
    incucyteDf = incucyteDf.astype(float).T

    return incucyteDf,plateLayoutDf
