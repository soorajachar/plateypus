#!/usr/bin/env python3 
import json,pickle,math,matplotlib,sys,os,string
import shutil
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
if os.name == 'nt':
    dirSep = '\\'
else:
    dirSep = '/'
from plateypus.dataprocessing.miscFunctions import cleanUpFlowjoCSV,reorderDfByInputOrder
from plateypus.dataprocessing import cytokineDataProcessing,singleCellDataProcessing,cellDataProcessing

dataTypeLevelNames = {'cyt':['Cytokine'],'cell':['CellType','Marker','Statistic'],'prolif':['Statistic'],'singlecell':['CellType']}
dataTypeDataFrameFileNames = {'cyt':'cytokineConcentrationPickleFile','cell':'cellStatisticPickleFile','prolif':'proliferationStatisticPickleFile','singlecell':'initialSingleCellPickleFile','killing':'killingIndexPickleFile'}
plateRowLetters = list(string.ascii_uppercase)[:16]
plateColumnNumbers = list(range(1,25))

def returnMultiIndex(sortedData,sortedFiles,dataType,folderName):
    if(dataType == 'cyt'):
        newMultiIndex = cytokineDataProcessing.parseCytokineCSVHeaders(pd.read_csv('inputData'+dirSep+'bulkCSVFiles'+dirSep+'A1_'+dataType+'.csv').columns)
    elif(dataType == 'cell'):
        if 'antibodyPanel-'+folderName+'.csv' in os.listdir('misc'):
            panelData = pd.read_csv('misc'+dirSep+'antibodyPanel-'+folderName+'.csv',)
            newMultiIndex = cellDataProcessing.parseCellCSVHeaders(pd.read_csv('inputData'+dirSep+'bulkCSVFiles'+dirSep+'A1_'+dataType+'.csv').columns,panelData=panelData)
        else:
            newMultiIndex = cellDataProcessing.parseCellCSVHeaders(pd.read_csv('inputData'+dirSep+'bulkCSVFiles'+dirSep+'A1_'+dataType+'.csv').columns)
    elif(dataType == 'cytcorr'):
        newMultiIndex = []
    if dataType != 'singlecell':
        return sortedData,newMultiIndex
    else:
        return sortedFiles,newMultiIndex

def decodeBarcodedPlates(experimentParameters,folderName,dataType,reverse):
    path = 'inputData'+dirSep+'bulkCSVFiles'+dirSep
    barcodingDict = experimentParameters['barcodingDict']
    bIndex  = 0
    alignmentColumnOrder = []
    for barcodedPlate in barcodingDict:
        barcodedCSV = pd.read_csv(path+barcodedPlate+'_'+dataType+'.csv')
        for decodedPlate in barcodingDict[barcodedPlate]:
            decodingIndices = [0]
            barcodes = barcodingDict[barcodedPlate][decodedPlate]
            for col,columnHeader in enumerate(barcodedCSV.columns):
                includeInDecoding = True
                for barcode in barcodes:
                    if barcode.lower().replace(' ','') not in columnHeader.lower().replace(' ',''):
                        includeInDecoding = False
                if includeInDecoding:
                    decodingIndices.append(col)
            decodingIndices.append(-1)
            decodedCSV = barcodedCSV.iloc[:,decodingIndices]
            newColumns = []
            for col,column in enumerate(decodedCSV.columns):
                if col == 0:
                    newColumns.append('')
                elif col == len(decodedCSV.columns)-1:
                    newColumns.append('')
                else:
                    populationStatisticSplit = column.split(' | ')
                    population = populationStatisticSplit[0]
                    splitPopulations = population.split('/')
                    populationsToKeep = []
                    for i,splitPopulation in enumerate(splitPopulations):
                        keepPopulation = True
                        for barcode in barcodes:
                            if barcode.lower().replace(' ','') in splitPopulation.lower().replace(' ',''):
                                keepPopulation = False
                        if keepPopulation:
                            populationsToKeep.append(splitPopulation)
                    newColumn = ' | '.join(['/'.join(populationsToKeep),populationStatisticSplit[1]])
                    newColumns.append(newColumn)
            decodedCSV.columns = newColumns
            if bIndex == 0:
                alignmentColumnOrder = decodedCSV.columns.tolist()
            else:
                if alignmentColumnOrder != decodedCSV.columns.tolist():
                    #Realign columns
                    newDecodedCSV = decodedCSV.iloc[:,:-1][alignmentColumnOrder[:-1]]
                    decodedCSV.iloc[:,:-1] = newDecodedCSV.values
                    decodedCSV.columns = alignmentColumnOrder
            
            decodedCSV.to_csv(path+decodedPlate+'_'+dataType+'.csv',index=False)
            bIndex+=1
        
def unpackMultiplexedPlates(experimentParameters,folderName,dataType,reverse):
    #A1->A1,A2->A2,A3- >A1,B1->A4,B2->A3
    #1.1,1.3,1.5,3.1,3.3,3.5 -> A1; 1.2,1.4,1.6,3.2,3.4,3.6->A2; 2.1,2.3,3.5,4.1,4.3,4.5->A4; 2.2,2.4,2.6,4.2,4.4,4.6->A3
    #1-24
    #A-P
    wellIDConversionDict = {}
    for i,plateLetter in enumerate(plateRowLetters):
        for j,plateNumber in enumerate(plateColumnNumbers):
            wellIDConversionDict[plateLetter+str(plateNumber)] = plateRowLetters[int(i/2)]+str(plateColumnNumbers[int(j/2)])
    plateIDConversionDict = {}
    for multiplexedPlateName in experimentParameters['unpackingDict']:
        for i,multiplexedWellPos in enumerate(experimentParameters['unpackingDict'][multiplexedPlateName]):
            wellIDs = []
            if multiplexedWellPos != '':
                if i < 2:
                    j = 0
                else:
                    j = 1
                for rowIndex in range(j,len(plateRowLetters),2):
                    for colIndex in range(i%2,len(plateColumnNumbers),2):
                        wellID = plateRowLetters[rowIndex]+str(plateColumnNumbers[colIndex])
                        wellIDs.append(wellID)
            plateIDConversionDict[multiplexedWellPos] = wellIDs
    #Specimen_001_A1_A01_001.fcs
    multiplexedPlateNames = experimentParameters['unpackingDict'].keys()
    #Unpacking
    if not reverse:
        for multiplexedPlateName in multiplexedPlateNames:
            multiplexedWellPoses = experimentParameters['unpackingDict'][multiplexedPlateName]
            with open('inputData'+dirSep+'bulkCSVFiles'+dirSep+multiplexedPlateName+'_'+dataType+'.csv', 'r') as f:
                multiplexedCSVLines = f.readlines()
            for multiplexedWellPos in multiplexedWellPoses:
                if multiplexedWellPos != '':
                    wellIDsInThisPos = plateIDConversionDict[multiplexedWellPos]
                    linesToMove = []
                    newCSVLines = []
                    for lineNum,line in enumerate(multiplexedCSVLines):
                        if lineNum in [0,len(multiplexedCSVLines)-2,len(multiplexedCSVLines)-1]:
                            newCSVLines.append(line)
                        else:
                            #CyTEK Explorer
                            if ' Well' in line:
                                wellID = line.split(' ')[0]
                            #FACSDiva
                            else:
                                fileName = line.split(',')[0]
                                wellID = fileName.split('_')[2]
                            #Well ID is in pos or first line or last two lines
                            if wellID in wellIDsInThisPos:
                                #underscorePoses = [pos for pos, char in enumerate(line) if char == '_']
                                newWellID = wellIDConversionDict[wellID]
                                newLine = line.replace(wellID,newWellID)
                                newCSVLines.append(newLine)
                    with open('inputData'+dirSep+'bulkCSVFiles'+dirSep+multiplexedWellPos+'_'+dataType+'.csv', 'w') as f:
                        for item in newCSVLines:
                            f.write("%s" % item)
    #Repacking
    else:
        #Invert well conversion dictionary
        for mpos,multiplexedPlateName in enumerate(experimentParameters['unpackingDict']):
            newCSVLines = []
            for mpos2,multiplexedWellPos in enumerate(experimentParameters['unpackingDict'][multiplexedPlateName]):
                specificWellIDConversionDict = {}
                for well in plateIDConversionDict[multiplexedWellPos]:
                    specificWellIDConversionDict[well] = wellIDConversionDict[well]
                invertedWellIDConversionDict = {v: k for k, v in specificWellIDConversionDict.items()}
                if multiplexedWellPos != '':
                    #Read unpacked csv, rename to name-unpacked.csv, delete old name (as repacked csv names can overlap with these)
                    with open('inputData'+dirSep+'bulkCSVFiles'+dirSep+multiplexedWellPos+'_'+dataType+'.csv', 'r') as f:
                        multiplexedCSVLines = f.readlines()
                    shutil.copyfile('inputData'+dirSep+'bulkCSVFiles'+dirSep+multiplexedWellPos+'_'+dataType+'.csv','inputData'+dirSep+'bulkCSVFiles'+dirSep+multiplexedWellPos+'-unpacked_'+dataType+'.csv')
                    os.remove('inputData'+dirSep+'bulkCSVFiles'+dirSep+multiplexedWellPos+'_'+dataType+'.csv')
                    
                    for lineNum,line in enumerate(multiplexedCSVLines):
                        #If first plate to repack, use first line
                        #If last plate to repack, use last 2 lines
                        if lineNum in [0,len(multiplexedCSVLines)-2,len(multiplexedCSVLines)-1]:
                            if (lineNum == 0 and mpos2 == 0) or (lineNum in [len(multiplexedCSVLines)-2,len(multiplexedCSVLines)-1] and mpos2 == len(experimentParameters['unpackingDict'][multiplexedPlateName])-1):
                                newCSVLines.append(line)
                        else:
                            #CyTEK Explorer
                            if ' Well' in line:
                                wellID = line.split(' ')[0]
                            #FACSDiva
                            else:
                                fileName = line.split(',')[0]
                                wellID = fileName.split('_')[2]
                            newWellID = invertedWellIDConversionDict[wellID]        
                            #print(wellID+'->'+newWellID)
                            newLine = line.replace(wellID,newWellID)
                            newCSVLines.append(newLine)
                     
            with open('inputData'+dirSep+'bulkCSVFiles'+dirSep+multiplexedPlateName+'_'+dataType+'.csv', 'w') as f:
                for item in newCSVLines:
                    f.write("%s" % item)

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
    shutil.move('inputData'+dirSep+'bulkCSVFiles'+dirSep+'temp-'+fileName,'inputData'+dirSep+'bulkCSVFiles'+dirSep+fileName)

def createBaseDataFrame(experimentParameters,folderName,experimentNumber,dataType,layoutDict):
    if experimentParameters['format'] == 'tube':
        fullFormatDf = pickle.load(open('misc'+dirSep+'tubeLayout-'+folderName+'-cell.pkl','rb'))
        dfList = []
        levelLabelDict = experimentParameters['levelLabelDict']
        for fileName in os.listdir('inputData'+dirSep+'bulkCSVFiles'+dirSep):
            if '.csv' in fileName:
                performCommaCheck(fileName)
                
                bulkTubeCSVFileName = fileName
                columnMultiIndexTuples = cellDataProcessing.parseCellCSVHeaders(pd.read_csv('inputData'+dirSep+'bulkCSVFiles'+dirSep+bulkTubeCSVFileName).columns)
                columnMultiIndex = pd.MultiIndex.from_tuples(columnMultiIndexTuples,names=['CellType','Marker','Statistic'])

                fullData = pd.read_csv('inputData'+dirSep+'bulkCSVFiles'+dirSep+bulkTubeCSVFileName,header=0)
                if 'Unnamed' in fullData.columns[-1]:
                    data = fullData.iloc[:-2,1:-1].values
                else:
                    data = fullData.iloc[:-2,1:].values
                sampleNames = fullData.iloc[:-2,0].values.ravel()
                sampleIndexStart = fullFormatDf.values.ravel().tolist().index(sampleNames[0])
                sampleIndexEnd = fullFormatDf.values.ravel().tolist().index(sampleNames[-1])
                rowMultiIndex = fullFormatDf.iloc[sampleIndexStart:sampleIndexEnd+1,:].index
                
                timeDataList = []
                timeSubsets = []
                times = levelLabelDict[list(levelLabelDict.keys())[-1]]

                #Can use sample name file to assign time values
                if 'sampleNameFile.xlsx' in os.listdir('misc') or 'sampleNameFile.csv' in os.listdir('misc'):
                    if 'sampleNameFile.xlsx' in os.listdir('misc'): 
                        sampleNameDf = pd.read_excel('misc'+dirSep+'sampleNameFile.xlsx')
                    else:
                        sampleNameDf = pd.read_csv('misc'+dirSep+'sampleNameFile.csv')
                    if 'Time' in sampleNameDf.columns:
                        for time in times:
                            timeIndices = []
                            for row in range(sampleNameDf.shape[0]):
                                if sampleNameDf[list(levelLabelDict.keys())[-1]].values[row] == time:
                                    timeIndices.append(row)
                            timeSubsets.append(timeIndices)
                    #Otherwise just assume 1 timepoint (HACK NEED TO FIX EVENTUALLY)
                    else:
                        timeSubsets.append(list(range(data.shape[0])))
                #Otherwise just assume 1 timepoint (HACK NEED TO FIX EVENTUALLY)
                else:
                    timeSubsets.append(list(range(data.shape[0])))

                for timeSubset in timeSubsets: 
                    dataList = []
                    columnTupleList = []
                    for i,columnTuple in enumerate(columnMultiIndexTuples):
                        ser = pd.Series(data[timeSubset,i],index=rowMultiIndex)
                        dataList.append(ser)
                        columnTupleList.append(tuple(columnTuple))
                    fullExperimentDf = pd.concat(dataList,keys=columnTupleList,names=['CellType','Marker','Statistic'])
                    timeDataList.append(fullExperimentDf)
                
                k = pd.concat(timeDataList,keys=times,names=[list(levelLabelDict.keys())[-1]])
                repeatList = []
                for name in k.index:
                    if name not in repeatList:
                        repeatList.append(name)
                    else:
                        print('Repeated:')
                        print(name)
                partialExperimentDf = pd.concat(timeDataList,keys=times,names=[list(levelLabelDict.keys())[-1]]).unstack(list(levelLabelDict.keys())[-1])
                dfList.append(partialExperimentDf)

        fullExperimentDf = pd.concat(dfList)
    else:
        realDataType = dataType
        
        #Reverse order; first decode 96 well barcoded plates, then repack the unbarcoded 96 well plates to a 384 well plate 
        if experimentParameters['overallPlateDimensions'][0] == 16 and 'unpackingDict' in list(experimentParameters.keys()):
            reverse = True
        #Normal order; first decode 384 well barcoded plates, then unpack the unbarcoded 384 well plates to 96 well plates
        else:
            reverse = False
        if 'barcodingDict' in list(experimentParameters.keys()):
            decodeBarcodedPlates(experimentParameters,folderName,dataType,reverse)
        if 'unpackingDict' in list(experimentParameters.keys()):
            unpackMultiplexedPlates(experimentParameters,folderName,dataType,reverse)

        #Legacy experiment parameter files compatibility
        if 'paired' in experimentParameters.keys():
            if experimentParameters['paired']:
                numRowPlates = 2
            else:
                numRowPlates = 1
            numColumnPlates = int(experimentParameters['numPlates']/numRowPlates)
        else:
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
        plateDimensions = experimentParameters['overallPlateDimensions'] 
        levels = list(levelLabelDict.keys())
        conditionLevels = list(levelLabelDict.keys())[:-1]
        conditionLevelValues = levelLabelDict.copy()
        del conditionLevelValues[list(levelLabelDict.keys())[-1]]
        allLevelValues = experimentParameters['levelLabelDict']
        
        sortedData,sortedFiles = cleanUpFlowjoCSV(plateNames,folderName,dataType,experimentParameters)
        allRawData,newLevelList = returnMultiIndex(sortedData,sortedFiles,realDataType,folderName)

        #print(allRawData)
        dfList = []
        for rawData,plateID in zip(allRawData,plateNames):
            fullTupleList = []
            index = 0
            for row in range(rawData.shape[0]):
                sampleID = int(rawData.iloc[row,0])-1
                wellID = plateRowLetters[int(sampleID/plateDimensions[1])]+str(plateColumnNumbers[sampleID % plateDimensions[1]])
                fullID = plateID+'-'+wellID
                
                sampleLocation = np.argwhere(identificationMatrix == fullID)[0]
                column = []
                tupleList = []
                for levelID in layoutDict['keys']:
                    levelValueID = layoutDict['keys'][levelID][sampleLocation[0],sampleLocation[1]]
                    if levelValueID != 'blank':
                        level = levels[levelID]
                        levelValue = allLevelValues[level][levelValueID]
                        tupleList.append(levelValue)
                if len(tupleList) != 0:
                    fullTupleList.append(tupleList)
                index+=1
            mi = pd.MultiIndex.from_tuples(fullTupleList,names=levels)
            columnSeriesList = []
            columnTupleList = []
            for column,columnTuple in enumerate(newLevelList):
                columnSeries = pd.Series(rawData.values[:,column+1],index=mi)
                columnSeriesList.append(columnSeries)
                columnTupleList.append(tuple(columnTuple))
            #print(columnTupleList)
            #print(columnSeriesList)
            #print(rawData)
            #print(newLevelList)
            platedf = pd.concat(columnSeriesList,axis=0,keys=columnTupleList,names=dataTypeLevelNames[realDataType])
            dfList.append(platedf)

        idx=pd.IndexSlice 
        fullExperimentDf = pd.concat(dfList)
        #dfl = [fullExperimentDf.xs([12.0],level=['Time']),fullExperimentDf.xs([60.0],level=['Time']),fullExperimentDf.xs([96.0],level=['Time']),fullExperimentDf.xs([156.0],level=['Time'])]
        tempdf = fullExperimentDf.to_frame('temp')
        temp = []
        for row in range(fullExperimentDf.shape[0]):
            name = list(tempdf.iloc[row,:].name)
            if name in temp:
                print(name)
                print(row)
            else:
                temp.append(name)
        #Remove blanks
        for i,level in enumerate(fullExperimentDf.index.names):
            tempLevelValues = pd.unique(fullExperimentDf.index.get_level_values(level))
            if 'Blank' in tempLevelValues.tolist():
                fullExperimentDf = fullExperimentDf.drop('Blank',level=i)
        temp = []
        temp2 = []
        tempdf = fullExperimentDf.to_frame('wat')
        for row in range(fullExperimentDf.shape[0]):
            name = list(tempdf.iloc[row,:].name)
            if name not in temp:
                temp.append(name)
            else:
                temp2.append(name)
        fullExperimentDf = fullExperimentDf.unstack(list(levelLabelDict.keys())[-1])
    
    fullExperimentDf = reorderDfByInputOrder(experimentParameters,fullExperimentDf)
    return fullExperimentDf

def convertDataFramesToExcel(folderName,secondPath,dataType,df):
    writer = pd.ExcelWriter('outputData'+dirSep+'excelFiles'+dirSep+'excelFile-'+folderName+'-'+dataType+'.xlsx')
    if dataType == 'cyt':
        dfg = pickle.load(open('outputData'+dirSep+'pickleFiles'+dirSep+'cytokineGFIPickleFile-'+folderName+'.pkl','rb'))
        dfc = pickle.load(open('outputData'+dirSep+'pickleFiles'+dirSep+dataTypeDataFrameFileNames[dataType]+'-'+folderName+'.pkl','rb'))
        dfg.to_excel(writer,'MFI')
        dfc.to_excel(writer,'Concentration')
    elif dataType == 'killing':
        dfk = pickle.load(open('outputData'+dirSep+'pickleFiles'+dirSep+dataTypeDataFrameFileNames[dataType]+'-'+folderName+'.pkl','rb'))
        dfk.to_excel(writer,'Killing Index')
    else:
        for statistic in list(pd.unique(df.index.get_level_values('Statistic'))):
            statisticDf = df.xs(statistic,level='Statistic')
            statisticDf.to_excel(writer,statistic)
    writer.save()

def saveFinalDataFrames(folderName,secondPath,experimentNumber,dataType,fullExperimentDf,excel_data):
    fullExperimentDf = fullExperimentDf.astype(float)
    with open('outputData'+dirSep+'pickleFiles'+dirSep+dataTypeDataFrameFileNames[dataType]+'-'+folderName+'.pkl', "wb") as f:
        pickle.dump(fullExperimentDf, f)
    convertDataFramesToExcel(folderName,secondPath,dataType,fullExperimentDf)
    print(fullExperimentDf)
