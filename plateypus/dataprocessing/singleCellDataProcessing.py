#! /usr/bin/env python3
#!/usr/bin/env python3 
import pickle,os,sys,string,shutil
import numpy as np
import pandas as pd
import tkinter as tk
import os
if os.name == 'nt':
    dirSep = '\\'
else:
    dirSep = '/'
from plateypus.dataprocessing.miscFunctions import reindexDataFrame,printProgressBar,extractValues,reorderDfByInputOrder,returnSpecificExtensionFiles
idx = pd.IndexSlice

plateRowLetters = list(string.ascii_uppercase)[:16]
plateColumnNumbers = list(range(1,25))

def produceSingleCellHeaders(cellTypes):
    newMultiIndexList = []
    for cellType in cellTypes:
        newMultiIndexList.append([cellType])
    return newMultiIndexList

def grabCellTypeList(experimentParameters):
    path = 'inputData'+dirSep+'singleCellCSVFiles'+dirSep
    if experimentParameters['format'] == 'plate':
        nrows = experimentParameters['overallPlateDimensions'][0]
        ncols = experimentParameters['overallPlateDimensions'][1]
        if 'unpackingDict' in experimentParameters.keys():
            orderingList = [str(x).zfill(3) for x in range(1,385)]
        else:
            orderingList = [str(x).zfill(3) for x in range(1,nrows*ncols+1)]

        #Plate case
        #Walk through a plate, grab each unique cell type
        cellTypeList = []
        folder = ''
        for fileName in os.listdir(path):
            if '.DS' not in fileName:
                if 'barcodingDict' in experimentParameters:
                    if fileName not in experimentParameters['barcodingDict']:
                        folder = fileName
                        break
                else:
                    folder = fileName
                    break

        for fileName in os.listdir(path+folder+dirSep):
            if '.DS' not in fileName:
                splitFile = fileName.split('_')
                parsingVal = ''
                for split in splitFile:
                    if split in orderingList:
                        parsingVal = split
                parsingPose = fileName.rindex(parsingVal)+4
                cellType = fileName[parsingPose:].split('.')[0].replace('/','\/')
                if cellType not in cellTypeList:
                    cellTypeList.append(cellType)
    else:
        cellTypeList = []
        sampleNameDf = pd.read_pickle('misc'+dirSep+'tubeLayout-'+os.getcwd().split(dirSep)[-1]+'-cell.pkl')
        fullSampleFileName = sampleNameDf.iloc[0,0]
        dotIndex = fullSampleFileName.rfind('.')
        sampleFileName = fullSampleFileName[:dotIndex]
        for fileName in os.listdir(path):
            if '.DS' not in fileName and sampleFileName in fileName:
                splitFile = fileName.split(sampleFileName[1:]+'_')
                cellType = splitFile[1].split('.')[0]
                if cellType not in cellTypeList:
                    cellTypeList.append(cellType)

    return cellTypeList

def createTubeSingleCellDataFrame(folderName,experimentParameters,fileNameDf):
    
    path = 'inputData'+dirSep+'singleCellCSVFiles'+dirSep
    
    fileNameDf = fileNameDf.stack().to_frame('fileName')

    #Grab common prefix for files
    fullFileName = ''
    for fileName in os.listdir(path):
        if '.DS' not in fileName:
            fullFileName = fileName
            break
    prefix = fullFileName.split('_')[0]

    #Remove extraneous decorations in marker names
    newColumns = []
    fcsDf = pd.read_csv(path+fullFileName,header=0)
    for column in fcsDf.columns:
        if '::' in column:
            column = column.split(' :: ')[1]
        #CyTOF
        if column[0].isdigit():
            column = column.split('_')[1]
        newColumns.append(column)

    cellTypeList = grabCellTypeList(experimentParameters)
    
    completeDfList = []
    for cellType in cellTypeList:
        for row in range(fileNameDf.shape[0]):
            fullFileName = fileNameDf.iloc[row,0]
            dotIndex = fullFileName.rfind('.')
            fileName = fullFileName[:dotIndex]
            trueFileName = prefix+'_'+fileName+'_'+cellType+'.csv'
            levelValues = list(fileNameDf.iloc[row,:].name) 
            fcsDf = pd.read_csv(path+trueFileName,header=0)
            eventNumber = fcsDf.shape[0]
            eventList = range(1,eventNumber+1)
            allLevelValues = []
            for event in eventList:
                allLevelValues.append([cellType]+levelValues+[event])
            newMultiIndex = pd.MultiIndex.from_tuples(allLevelValues,names=['CellType']+list(fileNameDf.index.names)+['Event'])
            newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=newColumns)
            completeDfList.append(newDf)
            printProgressBar(row + 1, fileNameDf.shape[0], prefix = ' Concatenating samples:', suffix = 'Complete', length = 50)

    completeDataFrame = pd.concat(completeDfList)
    completeDataFrame.columns.name = 'Marker'

    #Remove extraneous markers (namely -h parameters)
    columnsToKeep = []
    for col,column in enumerate(completeDataFrame.columns):
        if '-H' not in column and 'Time' not in column and '-W' not in column and column != 'GFP' and column != 'LiveDead':
            columnsToKeep.append(col)
    completeDataFrame = completeDataFrame.iloc[:,columnsToKeep]

    completeDataFrame.to_hdf('outputData'+dirSep+'pickleFiles'+dirSep+'initialSingleCellDf-channel-'+folderName+'.h5', key='df', mode='w')
    print(completeDataFrame)

def debarcodeSingleCellData(experimentParameters):
    barcodingDict = experimentParameters['barcodingDict']
    for barcodedPlateName in barcodingDict:
        inputPath = 'inputData'+dirSep+'singleCellCSVFiles'+dirSep+barcodedPlateName+dirSep
        for debarcodedPlateName in barcodingDict[barcodedPlateName]:
            #Create folder if not there; if there, delete folder before proceeding (avoids overlap)
            outputPath = 'inputData'+dirSep+'singleCellCSVFiles'+dirSep+debarcodedPlateName+dirSep
            if debarcodedPlateName in os.listdir('inputData'+dirSep+'singleCellCSVFiles'+dirSep):
                shutil.rmtree(outputPath[:-1])
            os.mkdir(outputPath[:-1])
            
            barcodes = barcodingDict[barcodedPlateName][debarcodedPlateName]
            #Search through all files in barcoded plate folder, move to each debarcoded plate folder
            for k,fileName in enumerate([x for x in os.listdir(inputPath) if '.DS' not in x]):
                allBarcodesInFileName = True
                for barcode in barcodes:
                    if barcode not in fileName:
                        allBarcodesInFileName = False
                        break
                if allBarcodesInFileName:
                    #Copy file to debarcoded plate folder
                    shutil.copyfile(inputPath+fileName,outputPath+fileName)
                    #Rename file; if __ exists, use suffix after __ as cell population name. otherwise no suffix
                    #In either case, remove barcode suffix
                    if '__' in fileName:
                        population = fileName.split('.')[0].split('__')[1]
                        newFileName = '_'.join(fileName.split('.')[0].split('__')[0].split('_')[:6]+[population])+'.csv'
                    else:
                        newFileName = '_'.join(fileName.split('.')[0].split('_')[:6]+['allCells'])+'.csv'
                    shutil.move(outputPath+fileName,outputPath+newFileName)
                printProgressBar(k + 1, len([x for x in os.listdir(inputPath) if '.DS' not in x]), prefix = ' Debarcoding '+barcodedPlateName+' to '+debarcodedPlateName+':', suffix = 'Complete', length = 50)

#If multiplexing option chosen
def demultiplexSingleCellData(experimentParameters):
    unpackingDict = experimentParameters['unpackingDict']
    unpackingPositionDict = {(0,0):0,(0,1):1,(1,0):2,(1,1):3}
    cellTypeList = grabCellTypeList(experimentParameters)

    #Currently in format A1; A2, B1; B2
    #Need to change to format A1: B1, A2, B2
    plateRowLetters = string.ascii_uppercase[:16]
    plateColumnNumbers = list(range(1,25))

    wellPlateRowLetters = string.ascii_uppercase[:8]
    wellPlateColumnNumbers = list(range(1,13))

    #Unpacking
    if experimentParameters['overallPlateDimensions'][0] == 8:
        #"multiplexingOption": "96->384 well", "unpackingDict": {"A1-2_B1-2": ["A1", "A2", "B2", "B1"]
        #Create appropriate folders in each population
        for combinedPlateName in list(unpackingDict.keys()):
            for unpackedPlateName in unpackingDict[combinedPlateName]:
                if unpackedPlateName != '':
                    if unpackedPlateName not in returnSpecificExtensionFiles('inputData'+dirSep+'singleCellCSVFiles'+dirSep,'',False):
                        os.mkdir('inputData'+dirSep+'singleCellCSVFiles'+dirSep+unpackedPlateName)
        fileNameDict = {}
        combinedPlateNames = list(unpackingDict.keys())
        for combinedPlateName in combinedPlateNames:
            #scale_Specimen_001_P9_P09_369_TCells.csv
            allFileNames = returnSpecificExtensionFiles('inputData'+dirSep+'singleCellCSVFiles'+dirSep+combinedPlateName,'',False)
            unpackedPlateNames = unpackingDict[combinedPlateName]
            for k,fileName in enumerate(allFileNames):
                sampleID = fileName.split('_')[3]
                currentRowLetter = sampleID[0]
                currentColumnNumber = int(sampleID[1:])
                #Get index of current row and column position
                currentRowLetterIndex = plateRowLetters.index(currentRowLetter)
                currentColumnNumberIndex = plateColumnNumbers.index(currentColumnNumber)
                #Demultiplex sample ids 
                wellPlateRowLetter = wellPlateRowLetters[int(currentRowLetterIndex/2)]
                wellPlateColumnNumber = wellPlateColumnNumbers[int(currentColumnNumberIndex/2)]
                newSampleID = str(wellPlateRowLetter)+str(wellPlateColumnNumber)
                newSampleID2 = str(wellPlateRowLetter)+str(wellPlateColumnNumber).zfill(2)
                newFileName = '_'.join(['_'.join(fileName.split('_')[:3]),newSampleID,newSampleID2,'_'.join(fileName.split('_')[-2:])])
                unpackedPlateIndex = unpackingPositionDict[(currentRowLetterIndex%2,currentColumnNumberIndex%2)]
                unpackedFolder = unpackedPlateNames[unpackedPlateIndex]
                trueFileName = newFileName
                #print(fileName+'->'+newFileName)
                fileNameDict['_'.join(trueFileName.split('_')[1:-1])] = '_'.join(fileName.split('_')[1:-1])
                completeNewFileName = unpackedFolder+dirSep+trueFileName
                shutil.copyfile('inputData'+dirSep+'singleCellCSVFiles'+dirSep+combinedPlateName+dirSep+fileName,'inputData'+dirSep+'singleCellCSVFiles'+dirSep+completeNewFileName)
                printProgressBar(k + 1, len(allFileNames), prefix = ' Unpacking '+combinedPlateName+':', suffix = 'Complete', length = 50)
    #Repacking
    else:
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
        fileNameDict = {}
        combinedPlateNames = list(unpackingDict.keys())
        for combinedPlateName in combinedPlateNames:
            if combinedPlateName in os.listdir('inputData'+dirSep+'singleCellCSVFiles'+dirSep):
                suffix = '-unpacked'
                shutil.move('inputData'+dirSep+'singleCellCSVFiles'+dirSep+combinedPlateName,'inputData'+dirSep+'singleCellCSVFiles'+dirSep+combinedPlateName+'-unpacked')
            else:
                suffix = ''
            os.mkdir('inputData'+dirSep+'singleCellCSVFiles'+dirSep+combinedPlateName)
            for plateName in unpackingDict[combinedPlateName]:
                if plateName == combinedPlateName:
                    suffix = '-unpacked'
                else:
                    suffix = ''
                specificWellIDConversionDict = {}
                for well in plateIDConversionDict[plateName]:
                    specificWellIDConversionDict[well] = wellIDConversionDict[well]
                invertedWellIDConversionDict = {v: k for k, v in specificWellIDConversionDict.items()}
                #scale_Specimen_001_P9_P09_369_TCells.csv
                allFileNames = returnSpecificExtensionFiles('inputData'+dirSep+'singleCellCSVFiles'+dirSep+plateName+suffix,'',False)
                unpackedPlateNames = unpackingDict[combinedPlateName]
                for k,fileName in enumerate(allFileNames):
                    sampleID = fileName.split('_')[3]
                    newSampleID =  invertedWellIDConversionDict[sampleID]
                    newSampleID2 = newSampleID[0]+newSampleID[1:].zfill(3)
                    newFileName = '_'.join(['_'.join(fileName.split('_')[:3]),newSampleID,newSampleID2,'_'.join(fileName.split('_')[-2:])])
                    fileNameDict['_'.join(newFileName.split('_')[1:-1])] = '_'.join(fileName.split('_')[1:-1])
                    shutil.copyfile('inputData'+dirSep+'singleCellCSVFiles'+dirSep+plateName+suffix+dirSep+fileName,'inputData'+dirSep+'singleCellCSVFiles'+dirSep+combinedPlateName+dirSep+newFileName)
                    printProgressBar(k + 1, len(allFileNames), prefix = ' Repacking '+plateName+':', suffix = 'Complete', length = 50)

    with open('misc'+dirSep+'fileNameDict.pkl','wb') as f:
        pickle.dump(fileNameDict,f)

def createPlateSingleCellDataFrame(folderName,experimentParameters,levelLayout,useBlankWells):
    
    path = 'inputData'+dirSep+'singleCellCSVFiles'+dirSep
    #Change CyTEK file names to match BD file names
    for folder in os.listdir(path):
        if '.DS' not in folder:
            for fileName in os.listdir(path+folder):
                if '.DS' not in fileName:
                    if ' Well' in fileName:
                        newFileNameComponents = fileName.split(' Well')
                        newFileName = newFileNameComponents[0].split('_')[0]+'_Specimen_001_'+newFileNameComponents[0].split('_')[1]+'_'+newFileNameComponents[0].split('_')[1][0]+newFileNameComponents[0].split('_')[1][1:].zfill(3)+newFileNameComponents[1]
                        #print('singleCellCSVFiles/'+folder+'/'+fileName+'->'+'singleCellCSVFiles/'+folder+'/'+newFileName)
                        shutil.move(path+folder+dirSep+fileName,path+folder+dirSep+newFileName)
    
    if 'barcodingDict' in experimentParameters:
        debarcodeSingleCellData(experimentParameters)
    if 'unpackingDict' in experimentParameters:
        demultiplexSingleCellData(experimentParameters)
    
    nrows = experimentParameters['overallPlateDimensions'][0]
    ncols = experimentParameters['overallPlateDimensions'][1]

    if 'unpackingDict' in experimentParameters.keys():
        orderingList = [str(x).zfill(3) for x in range(1,385)]
    else:
        orderingList = [str(x).zfill(3) for x in range(1,nrows*ncols+1)]

    completeKeyMatrix = np.dstack(list(levelLayout['keys'].values()))
    unraveledKeyMatrix = np.reshape(completeKeyMatrix,(completeKeyMatrix.shape[0]*completeKeyMatrix.shape[1],completeKeyMatrix.shape[2]))
    unraveledBlankMatrix = levelLayout['blank'].ravel()
    #print(unraveledBlankMatrix)

    sampleIndex = pd.MultiIndex.from_arrays([levelLayout['plateID'].ravel(),levelLayout['wellID'].ravel()],names=['Plate','Well'])
    sampleKeyDf = pd.DataFrame(unraveledKeyMatrix,index=sampleIndex,columns=list(experimentParameters['levelLabelDict'].keys()))
    sampleDf = sampleKeyDf.copy()
    rowsToKeep = []
    
    #print(sampleKeyDf.iloc[:,3].values.ravel())
    for row in range(sampleDf.shape[0]):
        for col in range(sampleDf.shape[1]):
            level = list(experimentParameters['levelLabelDict'].keys())[col]
            levelValueIndex = sampleKeyDf.iloc[row,col]
            if unraveledBlankMatrix[row] == -1:
                levelValue = experimentParameters['levelLabelDict'][level][levelValueIndex]
                sampleDf.iloc[row,col] = levelValue
            else:
                sampleDf.iloc[row,col] = 'Blank'
    #Drop blanks
    sampleDf = sampleDf.query("Time != 'Blank'")
    sampleDf.to_excel('outputData'+dirSep+'excelFiles'+dirSep+'fcsLabelingKey.xlsx')
    
    cellTypeList = grabCellTypeList(experimentParameters)

    newSampleDfList = []
    for cellType in cellTypeList:
        #Create dict that relates sample file locations to plate/wellIDs
        orderValDict = {}
        for plateName in list(np.unique(levelLayout['plateID'])):
            orderValDict2 = {}
            if 'DS' not in plateName:
                for fileName in os.listdir(path+plateName+dirSep):
                    if 'DS' not in fileName:
                        splitFile = fileName.split('_')
                        parsingVal = ''
                        for split in splitFile:
                            if split in orderingList:
                                parsingVal = split
                        parsingPose = fileName.rindex(parsingVal)
                        currentCellType = fileName[parsingPose+4:].split('.')[0]
                        if currentCellType == cellType:
                            orderVal = fileName[parsingPose:parsingPose+3]
                            wellID = fileName[:parsingPose-1].split('_')[-2]
                            orderValDict2[wellID] = path+plateName+dirSep+fileName
            orderValDict[plateName] = orderValDict2
        #print(orderValDict)
        #sys.exit(0)
        sampleTupleList,sampleList = [],[]
        for row in range(sampleDf.shape[0]):
            sampleID = list(sampleDf.iloc[row,:].name)
            plateID = sampleID[0]
            wellID = sampleID[1]
            sampleFileName = orderValDict[plateID][wellID]
            sampleList.append(sampleFileName)
            if type(cellType) != list:
                if type(cellType) == tuple:
                    cellType = list(cellType)
                else:
                    cellType = [cellType]
            sampleTuple = cellType+sampleDf.loc[idx[plateID,wellID],:].values.tolist()
            sampleTupleList.append(sampleTuple)
        sampleMI = pd.MultiIndex.from_tuples(sampleTupleList,names=['CellType']+list(sampleDf.columns))
        newSampleDf = pd.DataFrame(sampleList,index=sampleMI,columns=['fileName'])
        newSampleDfList.append(newSampleDf)
    
    fileNameDf = pd.concat(newSampleDfList)

    bothBool = False
    fullFileName = fileNameDf.iloc[0,0]
    fcsDf = pd.read_csv(fullFileName,header=0)
    for column in fcsDf.columns:
        if '::' in column:
            bothBool = True
    
    completeDfList = []
    for row in range(fileNameDf.shape[0]):
        fullFileName = fileNameDf.iloc[row,0]
        levelValues = list(fileNameDf.iloc[row,:].name) 
        fcsDf = pd.read_csv(fullFileName,header=0)
        eventNumber = fcsDf.shape[0]
        if eventNumber == 0:
            if not useBlankWells:
                tk.messagebox.showerror("Error", "Filename:\n"+fullFileName+"\nhas no events. Please re-export this file and try again.")
                sys.exit(0)
            else:
                newMatrix = np.zeros([1,fcsDf.shape[1]])
                fcsDf = pd.DataFrame(newMatrix,columns=fcsDf.columns)
                eventNumber+=1
        eventList = range(1,eventNumber+1)
        allLevelValues = []
        for event in eventList:
            allLevelValues.append(levelValues+[event])
        newMultiIndex = pd.MultiIndex.from_tuples(allLevelValues,names=list(fileNameDf.index.names)+['Event'])
        if bothBool:
            newColumns = [x.split(' :: ')[1] if '::' in x else x for x in fcsDf.columns]
        else:
            newColumns = fcsDf.columns
        newDf = pd.DataFrame(fcsDf.values,index=newMultiIndex,columns=newColumns)
        completeDfList.append(newDf)
        printProgressBar(row + 1, fileNameDf.shape[0], prefix = ' Concatenating samples:', suffix = 'Complete', length = 50)
    
    #Remove all unnecessary folders (only keep input/output)
    foldersToKeep = []
    for folder in [x for x in os.listdir('inputData'+dirSep+'singleCellCSVFiles'+dirSep) if '.DS' not in x]:
        #Keep all folders
        if 'unpackingDict' in experimentParameters and 'barcodingDict' in experimentParameters and experimentParameters['overallPlateDimensions'][0] == 16:
            if folder in experimentParameters['barcodingDict'] or folder in experimentParameters['unpackingDict']:
                foldersToKeep.append(folder)
        else:
            foldersToKeep.append(folder)

    for folder in [x for x in os.listdir('inputData'+dirSep+'singleCellCSVFiles'+dirSep) if '.DS' not in x]:
        if folder not in foldersToKeep:
            shutil.rmtree('inputData'+dirSep+'singleCellCSVfiles'+dirSep+folder)

    completeDataFrame = pd.concat(completeDfList)
    completeDataFrame.columns.name = 'Marker'

    #Remove extraneous markers (namely -h/-w parameters)
    columnsToKeep = []
    for col,column in enumerate(completeDataFrame.columns):
        #if '-H' not in column and 'Time' not in column and '-W' not in column:
        if '-H' not in column and 'Time' not in column and '-W' not in column and column != 'GFP' and 'Live' not in column:
            columnsToKeep.append(col)
    completeDataFrame = completeDataFrame.iloc[:,columnsToKeep]

    completeDataFrame.to_hdf('outputData'+dirSep+'pickleFiles'+dirSep+'initialSingleCellDf-channel-'+folderName+'.h5', key='df', mode='w')
    print(completeDataFrame)
