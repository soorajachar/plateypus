#!/usr/bin/env python3
import pickle,os,json,math,subprocess
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk
from plateypus.setup.createTubeLayout import TubeLayoutPage 
from plateypus.setup.createPlateLayout import BlankSelectionPage

experimentParameters = {}
dataTypeDict = {'Cytokine':'cyt','Bulk/Single cell':'cell','Killing':'killing'}

class ExperimentSetupStartPage(tk.Frame):
    def __init__(self, master,fName,bPage):
        global folderName,backPage
        folderName = fName
        backPage = bPage
        tk.Frame.__init__(self, master)
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        v2 = tk.StringVar(value='inpt')
        rb2a = tk.Radiobutton(mainWindow, text="Enter experiment labels ",padx = 20, variable=v2, value='inpt')
        rb2b = tk.Radiobutton(mainWindow,text="Create experiment layout ",padx = 20, variable=v2, value='pl')
        rb2a.grid(row=0,column=0,sticky=tk.W)
        rb2b.grid(row=1,column=0,sticky=tk.W)
        
        #v3 = tk.StringVar(value='both')
        #l3 = tk.Label(mainWindow, text="for:     ")
        #l3.grid(row=0,column=1,rowspan=2,sticky=tk.EW)
        #rb3a = tk.Radiobutton(mainWindow, text="Cells",padx = 20, variable=v3, value='cell')
        #rb3b = tk.Radiobutton(mainWindow,text="Cytokines",padx = 20, variable=v3, value='cyt')
        #rb3c = tk.Radiobutton(mainWindow,text="Cells and Cytokines",padx = 20, variable=v3, value='both')
        #rb3a.grid(row=0,column=2,sticky=tk.W)
        #rb3b.grid(row=1,column=2,sticky=tk.W)
        #rb3c.grid(row=2,column=2,sticky=tk.W)
        
        labelList,cbList,cbVarList = [],[],[]
        for i,dataT in enumerate(list(dataTypeDict.keys())):
            l = tk.Label(mainWindow,text=dataT).grid(row=i,column=3,sticky=tk.W)
            v = tk.BooleanVar(value=False)
            cb = tk.Checkbutton(mainWindow, variable=v)
            cb.grid(row=i,column=1,sticky=tk.W)
            cbList.append(cb)
            cbVarList.append(v)
        
        def experimentLayout():
            for dataType in dataTypeList:
                if 'experimentParameters-'+folderName+'-'+dataType+'.json' in os.listdir('misc'):
                    experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'-'+dataType+'.json','r'))
                    break
            if 'format' not in experimentParameters.keys():
                experimentParameters['format'] = 'plate'
            if experimentParameters['format'] == 'plate':
                plateDimensions = experimentParameters['overallPlateDimensions']
                levelLabelDict = experimentParameters['levelLabelDict']
                conditionLevelValues = levelLabelDict.copy()
                del levelLabelDict.copy()[list(levelLabelDict.keys())[-1]]
                levelValues = []
                for level in experimentParameters['levelLabelDict']:
                    levelValues.append(experimentParameters['levelLabelDict'][level])
                maxNumLevelValues = len(max(levelValues,key=len))
                levels = list(experimentParameters['levelLabelDict'].keys())
                if 'killing' in dataTypeList:
                    levels = levels[:-1]
                    levelValues = levelValues[:-1]
                master.switch_frame(BlankSelectionPage,folderName,levels,levelValues,maxNumLevelValues,experimentParameters['numPlates'],plateDimensions,dataTypeList,ExperimentSetupStartPage,bPage)
            #Tube mode
            else:
                fcsFiles = []
                
                if 'A1_cell.csv' not in os.listdir('inputData/bulkCSVFiles/'):
                    if len(os.listdir('inputData/fcsFiles')) > 1:
                        for fcsName in os.listdir('inputData/fcsFiles'):
                            if '.fcs' in fcsName:
                                fcsFiles.append(fcsName)
                        if len(fcsFiles) == 0:
                            fcsFiles = ['                 ']
                    elif 'sampleNameFile.xlsx' in os.listdir('misc'):
                        fcsFiles = list(pd.read_excel('misc/sampleNameFile.xlsx',header=0).iloc[:,0].values.ravel())
                else:
                    bulkStatFile = pd.read_csv('inputData/bulkCSVFiles/A1_cell.csv')
                    for row in range(bulkStatFile.shape[0]):
                        if bulkStatFile.iloc[row,0] not in ['Mean','SD']:
                            fcsFiles.append(bulkStatFile.iloc[row,0])

                master.switch_frame(TubeLayoutPage,folderName,experimentParameters['levelLabelDict'],len(fcsFiles),v3.get(),ExperimentSetupStartPage,bPage)
        
        def collectInput():
            global dataTypeList
            dataTypeList = [list(dataTypeDict.values())[i] for i,x in enumerate(cbVarList) if x.get()]
            if v2.get() == 'inpt':
                master.switch_frame(ExperimentFormatPage,folderName)
            elif v2.get() == 'pl':
                experimentLayout()

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInput()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class ExperimentFormatPage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        v1 = tk.StringVar()
        l1 = tk.Label(mainWindow,text='Collection Format: ')
        rb2a = tk.Radiobutton(mainWindow, text="Tube",padx = 20, variable=v1, value='tube')
        rb2b = tk.Radiobutton(mainWindow,text="Plate",padx = 20, variable=v1, value='plate')
        v1.set('tube')
        
        l1.grid(row=0,column=0)
        rb2a.grid(row=1,column=0)
        rb2b.grid(row=2,column=0)
        
        def collectInputs():
            if v1.get() == 'tube':
                experimentParameters['format'] = 'tube'
                master.switch_frame(TubeExperimentParameterPage,folderName)
            else:
                experimentParameters['format'] = 'plate'
                master.switch_frame(PlateExperimentParameterPage,folderName)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ExperimentSetupStartPage,folderName,backPage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class TubeExperimentParameterPage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        

        l4 = tk.Label(mainWindow, text="How many samples are there in this experiment?")
        e0 = tk.Entry(mainWindow)
        l4.grid(row=0,column=0)
        e0.grid(row=0,column=1)
        
        l5 = tk.Label(mainWindow, text="How many unique levels are there in this experiment (excluding Time)?")
        e1 = tk.Entry(mainWindow)
        l5.grid(row=2,column=0)
        e1.grid(row=2,column=1)
        
        if 'sampleNameFile.xlsx' in os.listdir('misc'):
            df = pd.read_excel('misc/sampleNameFile.xlsx')
            columns = []
            for column in df:
                if column not in ['','FileName','Time']:
                    columns.append(column)
            e0.insert(tk.END,str(df.shape[0]))
            e1.insert(tk.END,str(len(columns)))

        def collectInputs():
            experimentParameters['numSamples'] = int(e0.get())
            experimentParameters['numAllLevels'] = int(e1.get())+1
            master.switch_frame(allLevelNamePage,folderName)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ExperimentFormatPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class PlateExperimentParameterPage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        
        Multiplexing_Options = {'96->384 Well Plate','Barcoding'}

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l7 = tk.Label(mainWindow, text="Multiplexing Options:")
        multiplex = tk.StringVar(value='None')
        plateMultiplexVar = tk.BooleanVar() 
        plateMultiplexVar2 = tk.BooleanVar() 
        barcodeMultiplexVar = tk.BooleanVar() 
        plateMultiplexCb = tk.Checkbutton(mainWindow,text='96->384 well',variable=plateMultiplexVar,onvalue=True,offvalue=False)
        plateMultiplexCbRev = tk.Checkbutton(mainWindow,text='384->96 well',variable=plateMultiplexVar2,onvalue=True,offvalue=False)
        plateMultiplexCbRev.config(state=tk.DISABLED)
        plateMultiplexCb.config(state=tk.NORMAL)
        barcodeMultiplexCb = tk.Checkbutton(mainWindow,text='Barcoding',variable=barcodeMultiplexVar,onvalue=True,offvalue=False)
        barcodingNumberLabel = tk.Label(mainWindow,text='# plates per barcoded plate') 
        barcodingNumberEntry = tk.Entry(mainWindow,width=4)

        v = tk.IntVar()
        v2 = tk.IntVar()
        v3 = tk.IntVar()
        v4 = tk.IntVar()
        
        def disable():
            if v3.get() == 384:
                plateMultiplexCb.config(state=tk.DISABLED)
                plateMultiplexCbRev.config(state=tk.NORMAL)
            else:
                plateMultiplexCbRev.config(state=tk.DISABLED)
                plateMultiplexCb.config(state=tk.NORMAL)

        l3 = tk.Label(mainWindow, text="""What format were the samples collected in?:""")
        rb3a = tk.Radiobutton(mainWindow, text="96 well plate",padx = 20, variable=v3, value=96,command=lambda: disable())
        rb3b = tk.Radiobutton(mainWindow,text="384 well plate",padx = 20, variable=v3, value=384,command=lambda: disable())
        l3.grid(row=0,column=0)
        rb3a.grid(row=0,column=1)
        rb3b.grid(row=0,column=2)
        v3.set(96)

        l4 = tk.Label(mainWindow, text="How many plates were used in the experiment?")
        e0 = tk.Entry(mainWindow)
        l4.grid(row=1,column=0)
        e0.grid(row=1,column=1)

        l6 = tk.Label(mainWindow, text="Enter the number of condition levels (excluding Time): ")
        e2 = tk.Entry(mainWindow)
        l6.grid(row=3,column=0)
        e2.grid(row=3,column=1)
        
        l7.grid(row=4,column=0)
        plateMultiplexCb.grid(row=4,column=1)
        plateMultiplexCbRev.grid(row=5,column=1)
        barcodeMultiplexCb.grid(row=4,column=2)
        barcodingNumberLabel.grid(row=4,column=3)
        barcodingNumberEntry.grid(row=4,column=4)
         
        def collectInputs():
            experimentParameters['numPlates'] = int(e0.get())
            
            #Incucyte datasets have time encoded in them already
            if 'killing' in dataTypeList:
                experimentParameters['numAllLevels'] = int(e2.get())
            else:
                experimentParameters['numAllLevels'] = int(e2.get())+1
            if v3.get() == 384:
                experimentParameters['overallPlateDimensions'] = [16,24]
            else:
                experimentParameters['overallPlateDimensions'] = [8,12]
            
            if (plateMultiplexVar.get() or plateMultiplexVar2.get()) and barcodeMultiplexVar.get():
                multiplexingOption = '96->384 well + Barcoding'
                numberOfBarcodes = barcodingNumberEntry.get()
            else:
                if plateMultiplexVar.get() or plateMultiplexVar2.get():
                    multiplexingOption = '96->384 well'
                    numberOfBarcodes = 1
                elif barcodeMultiplexVar.get():
                    multiplexingOption = 'Barcoding'
                    numberOfBarcodes = barcodingNumberEntry.get()
                else:
                    multiplexingOption = 'None'
                    numberOfBarcodes = 1
            if ',' in barcodingNumberEntry.get():
                numberOfBarcodes = list(map(int,numberOfBarcodes.split(',')))
            else:
                numberOfBarcodes = [int(numberOfBarcodes)]
            global globalMultiplexingVar
            globalMultiplexingVar = multiplexingOption
            experimentParameters['multiplexingOption'] = multiplexingOption
            if multiplexingOption == '96->384 well':
                master.switch_frame(MultiplexingPage,multiplexingOption)
            elif multiplexingOption == '96->384 well + Barcoding':
                master.switch_frame(BarcodingPage,multiplexingOption,numberOfBarcodes)
            elif multiplexingOption == 'Barcoding':
                master.switch_frame(BarcodingPage,multiplexingOption,numberOfBarcodes)
            else:
                master.switch_frame(allLevelNamePage,folderName)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ExperimentFormatPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class BarcodingPage(tk.Frame):
    def __init__(self, master,mpo,nb):
        tk.Frame.__init__(self, master)
        multiplexingOption = mpo
        BarcodingWindow = tk.Frame(self)
        BarcodingWindow.pack()
        numPl = experimentParameters['numPlates'] 
        
        numberOfBarcodes = math.ceil(np.log2(max(nb)))
        entryBarcodeLabelList = []
        for barcode in range(numberOfBarcodes):
            e = tk.Entry(BarcodingWindow,width=8)
            e.grid(row=0,column=2+barcode)
            e.insert(tk.END, 'Barcode '+str(barcode+1))
            entryBarcodeLabelList.append(e)
            #tk.Label(BarcodingWindow,text='Barcode '+str(barcode+1)).grid(row=0,column=2+barcode)
        
        if len(nb) == 1:
            numberOfBarcodedPlates = nb[0]
            if multiplexingOption == '96->384 well + Barcoding':
                if experimentParameters['overallPlateDimensions'][0] == 8: 
                    maxNumBarcoded = math.ceil(experimentParameters['numPlates'] / 4 / numberOfBarcodedPlates)
                    totalPlatesPerBarcodedPlate = 4*numberOfBarcodedPlates
                    plateSkip = 4
                #384->96
                else: 
                    maxNumBarcoded = math.ceil(4*experimentParameters['numPlates'] / numberOfBarcodedPlates)
                    totalPlatesPerBarcodedPlate = numberOfBarcodedPlates
                    plateSkip = 1 
            else:
                maxNumBarcoded = math.ceil(experimentParameters['numPlates'] / numberOfBarcodedPlates)
                totalPlatesPerBarcodedPlate = numberOfBarcodedPlates
                plateSkip = 1
            allIndividualPlateEntryList = []
            combinedPlateBarcodingVarList,combinedEntryList = [],[]
            for combinedPlateNum in range(maxNumBarcoded):
                combinedEntry = tk.Entry(BarcodingWindow,width=8)
                combinedEntry.insert(tk.END, 'A'+str(combinedPlateNum*totalPlatesPerBarcodedPlate+1)+'-'+str((combinedPlateNum+1)*totalPlatesPerBarcodedPlate))
                plateStart = combinedPlateNum*totalPlatesPerBarcodedPlate+1
                plateEnd = (combinedPlateNum+1)*totalPlatesPerBarcodedPlate
                allPlateIndices = list(range(plateStart,plateEnd+2,plateSkip))
                combinedEntry.grid(row=combinedPlateNum*numberOfBarcodedPlates+1,column=0)
                plateBarcodingVarList,plateEntryList = [],[]
                for plateNum in range(numberOfBarcodedPlates):
                    if plateNum == numberOfBarcodedPlates-1:
                        padvar = 10 
                    else:
                        padvar = 0
                    plateEntry = tk.Entry(BarcodingWindow,width=8)
                    if plateSkip == 4:
                        plateEntry.insert(tk.END, 'A'+str(allPlateIndices[plateNum])+'-'+str(allPlateIndices[plateNum+1]-1))
                    else:
                        plateEntry.insert(tk.END, 'A'+str(allPlateIndices[plateNum]))
                    plateEntry.grid(row=combinedPlateNum*numberOfBarcodedPlates+plateNum+1,column=1,pady=(0,padvar))
                    barcodingCbList,barcodingVarList = [],[]
                    for barcodeNum in range(numberOfBarcodes):
                        barcodeVar = tk.BooleanVar()
                        barcodeCb = tk.Checkbutton(BarcodingWindow,text='',variable=barcodeVar,onvalue=True,offvalue=False)
                        barcodeCb.grid(row=combinedPlateNum*numberOfBarcodedPlates+plateNum+1,column=2+barcodeNum,pady=(0,padvar))
                        barcodingCbList.append(barcodeCb)
                        barcodingVarList.append(barcodeVar)
                    plateBarcodingVarList.append(barcodingVarList)
                    plateEntryList.append(plateEntry)
                combinedPlateBarcodingVarList.append(plateBarcodingVarList)
                combinedEntryList.append(combinedEntry)
                allIndividualPlateEntryList.append(plateEntryList)
        else:
            if multiplexingOption == '96->384 well + Barcoding':
                if experimentParameters['overallPlateDimensions'][0] == 8: 
                    plateSkip = 4
                #384->96
                else: 
                    plateSkip = 1 
            else:
                plateSkip = 1
            allIndividualPlateEntryList = []
            combinedPlateBarcodingVarList,combinedEntryList = [],[]
            for elem,numberOfBarcodedPlates in enumerate(nb):
                summedPlateNums = sum(nb[:elem+1])
                plateEnd = plateSkip*summedPlateNums
                if elem == 0:
                    plateStart = 1
                    plateRow = 1
                else:
                    plateStart = sum(nb[:elem])*plateSkip + 1
                    plateRow = sum(nb[:elem]) + 1
                combinedEntry = tk.Entry(BarcodingWindow,width=8)
                combinedEntry.insert(tk.END, 'A'+str(plateStart)+'-'+str(plateEnd))
                combinedEntry.grid(row=plateRow,column=0)
                allPlateIndices = list(range(plateStart,plateEnd+2,plateSkip))
                plateBarcodingVarList,plateEntryList = [],[]
                for plateNum in range(numberOfBarcodedPlates):
                    if plateNum == numberOfBarcodedPlates-1 and plateNum != 0:
                        padvar = 10 
                    else:
                        padvar = 0
                    plateEntry = tk.Entry(BarcodingWindow,width=8)
                    if plateSkip == 4:
                        plateEntry.insert(tk.END, 'A'+str(allPlateIndices[plateNum])+'-'+str(allPlateIndices[plateNum+1]-1))
                    else:
                        plateEntry.insert(tk.END, 'A'+str(allPlateIndices[plateNum]))
                    plateEntry.grid(row=plateRow+plateNum,column=1,pady=(0,padvar))
                    barcodingCbList,barcodingVarList = [],[]
                    for barcodeNum in range(numberOfBarcodes):
                        barcodeVar = tk.BooleanVar()
                        barcodeCb = tk.Checkbutton(BarcodingWindow,text='',variable=barcodeVar,onvalue=True,offvalue=False)
                        barcodeCb.grid(row=plateRow+plateNum,column=2+barcodeNum,pady=(0,padvar))
                        barcodingCbList.append(barcodeCb)
                        barcodingVarList.append(barcodeVar)
                    plateBarcodingVarList.append(barcodingVarList)
                    plateEntryList.append(plateEntry)
                combinedPlateBarcodingVarList.append(plateBarcodingVarList)
                combinedEntryList.append(combinedEntry)
                allIndividualPlateEntryList.append(plateEntryList)

        l1 = tk.Label(BarcodingWindow,text='Barcoded Plate Name:')
        l2 = tk.Label(BarcodingWindow,text='Plate Name:')
        l1.grid(row=0,column=0)
        l2.grid(row=0,column=1)

        def collectInputs():
            barcodingDict = {}
            for combinedIndex,combinedEntry in enumerate(combinedEntryList):
                if combinedEntry.get() != '':
                    plateEntryList = allIndividualPlateEntryList[combinedIndex]
                    plateBarcodingDict = {}
                    for plateIndex,plateEntry in enumerate(plateEntryList):
                        if plateEntry.get() != '':
                            barcodingStringList = []
                            for barcodeIndex,barcodeVar in enumerate(combinedPlateBarcodingVarList[combinedIndex][plateIndex]):
                                if barcodeVar.get():
                                    barcodeSymbol = '+'
                                else:
                                    barcodeSymbol = '-'
                                barcodingStringList.append(entryBarcodeLabelList[barcodeIndex].get()+barcodeSymbol)
                            uncombinedPlateName = plateEntry.get()
                            if '_cyt' in plateEntry.get():
                                uncombinedPlateName = uncombinedPlateName[:-4]
                            elif '_cell' in plateEntry.get():
                                uncombinedPlateName = uncombinedPlateName[:-5]
                            plateBarcodingDict[uncombinedPlateName] = barcodingStringList
                    combinedPlateName = combinedEntry.get()
                    if '_cyt' in combinedEntry.get():
                        combinedPlateName = combinedPlateName[:-4]
                    elif '_cell' in combinedEntry.get():
                        combinedPlateName = combinedPlateName[:-5]
                    barcodingDict[combinedPlateName] = plateBarcodingDict
            experimentParameters['barcodingDict'] = barcodingDict
            if multiplexingOption == '96->384 well + Barcoding':
                master.switch_frame(MultiplexingPage,multiplexingOption)
            else:
                master.switch_frame(allLevelNamePage,folderName)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(PlateExperimentParameterPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class MultiplexingPage(tk.Frame):
    def __init__(self, master,mpo):
        tk.Frame.__init__(self, master)
        multiplexingOption = mpo
        Well384ConversionWindow = tk.Frame(self)
        Well384ConversionWindow.pack()
        
        reverse = False
        if multiplexingOption == '96->384 well + Barcoding':
            if experimentParameters['overallPlateDimensions'][0] == 16: 
                reverse = True

        if not reverse:
            maxNum384Plates = math.ceil(experimentParameters['numPlates'] / 4)
        else:
            maxNum384Plates = math.ceil(experimentParameters['numPlates'])
        Unpacking384PlateNameList = []
        Unpacking384List = []
        l1 = tk.Label(Well384ConversionWindow,text='Combined Plate Names: ')
        l1.grid(row=0,column=0)
        l2 = tk.Label(Well384ConversionWindow,text='Separated Plate Names: ')
        l2.grid(row=1,column=0)
        ttk.Separator().place(x=0, y=27, relwidth=1)
        ttk.Separator().place(x=157, y=0, relheight=0.65)
        for plateNum in range(maxNum384Plates):
            e1 = tk.Entry(Well384ConversionWindow,width=7)
            if not reverse:
                e1.insert(tk.END, 'A'+str(plateNum*4+1)+'-'+str((plateNum+1)*4))
            else:
                e1.insert(tk.END, 'A'+str(plateNum+1))
            e1.grid(row=0,column=plateNum*2+1)
            Unpacking384PlateNameList.append(e1)
            UnpackingWellPosList = []
            wellIndexList = [0,1,3,2]
            wellIndex = 0
            ttk.Separator().place(x=317+plateNum*160, y=0, relheight=0.65)
            for wellUnpackingRow in range(2):
                for wellUnpackingCol in range(2):
                    e2 = tk.Entry(Well384ConversionWindow,width=7)
                    e2.insert(tk.END, 'A'+str((plateNum*4)+wellIndexList[wellIndex]+1))
                    if wellUnpackingCol == 1:
                        padvar = 5
                    else:
                        padvar = 0
                    e2.grid(row=wellUnpackingRow+1,column=wellUnpackingCol+1+plateNum*2,padx=padvar)
                    UnpackingWellPosList.append(e2)
                    wellIndex+=1
            Unpacking384List.append(UnpackingWellPosList)

        def collectInputs():
            unpackingDict = {}
            for unpackingPlateNameEntry,unpackingWellPosList in zip(Unpacking384PlateNameList,Unpacking384List):
                tempList = []
                for unpackingEntry in unpackingWellPosList:
                    tempList.append(unpackingEntry.get())
                unpackingDict[unpackingPlateNameEntry.get()] = tempList
            experimentParameters['unpackingDict'] = unpackingDict
            master.switch_frame(allLevelNamePage,folderName)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(PlateExperimentParameterPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)
        
class allLevelNamePage(tk.Frame):
    def __init__(self, master,folderName):
        tk.Frame.__init__(self, master)
        numAllLevels = experimentParameters['numAllLevels']
            
        entryList1 = []
        entryList2 = []
        
        numericCheckBoxes = []
        numericCheckBoxVars = []

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        lt1 = tk.Label(mainWindow, text="Level Name").grid(row=0,column=1)
        lt2 = tk.Label(mainWindow, text="Number of Level Values").grid(row=0,column=2)
        lt3 = tk.Label(mainWindow, text="Numeric?").grid(row=0,column=3)
        
        if 'sampleNameFile.xlsx' in os.listdir('misc'):
            df = pd.read_excel('misc/sampleNameFile.xlsx')
            levels = []
            levelValueNums = []
            for level in df.columns:
                if level not in ['','FileName']:
                    levels.append(level)
                    levelValueNums.append(len(pd.unique(df[level])))
            j=0
        
        for conditionLevelNumber in range(1,numAllLevels+1):
            l1 = tk.Label(mainWindow, text="Condition "+str(conditionLevelNumber))
            e1 = tk.Entry(mainWindow)
            v = tk.IntVar()
            if conditionLevelNumber == 1 and 'killing' not in dataTypeList:
                #e1 = tk.Entry(mainWindow,state='readonly')
                e1.insert(tk.END, 'Time')
                e2 = tk.Entry(mainWindow)
                if 'sampleNameFile.xlsx' in os.listdir('misc'):
                    if 'Time' in df.columns:
                        e2.insert(tk.END,str(len(pd.unique(df['Time']))))
                    else: 
                        e2.insert(tk.END, '1')
                else:
                    e2.insert(tk.END, '1')
                v.set(1)
            else:
                e1 = tk.Entry(mainWindow)
                e2 = tk.Entry(mainWindow)
                v.set(0)
                if 'sampleNameFile.xlsx' in os.listdir('misc'):
                    e1.insert(tk.END, levels[j])
                    e2.insert(tk.END, str(levelValueNums[j]))
                    j+=1
                                    
            
            cb1 = tk.Checkbutton(mainWindow, text="",variable=v,onvalue=1,offvalue=0)
             
            l1.grid(row=conditionLevelNumber+1,column=0)
            e1.grid(row=conditionLevelNumber+1,column=1)
            e2.grid(row=conditionLevelNumber+1,column=2)
            entryList1.append(e1)
            entryList2.append(e2)
            cb1.grid(row=conditionLevelNumber+1,column=3)
            numericCheckBoxes.append(cb1)
            numericCheckBoxVars.append(v)
        
        def collectInputs():
            conditionNames = []
            numConditionLevelValues = []
            tiledLevels = []
            experimentParameters['columnVariableName'] = 'Time' 
            numericlevels = []
            for allLevelNumber in range(numAllLevels):
                #Remove column variable from condition name list
                if entryList1[allLevelNumber].get() == 'Time':
                    experimentParameters['numColumnLevelValues'] = int(entryList2[allLevelNumber].get())
                else:
                    conditionNames.append(str(entryList1[allLevelNumber].get()))
                    numConditionLevelValues.append(int(entryList2[allLevelNumber].get()))
                numericlevels.append(numericCheckBoxVars[allLevelNumber].get() == 1)
            if 'killing' in dataTypeList:
                experimentParameters['numColumnLevelValues'] = 0 
                experimentParameters['numConditionLevels'] = numAllLevels
            else:
                experimentParameters['numConditionLevels'] = numAllLevels - 1

            experimentParameters['conditionLevelNames'] = conditionNames
            experimentParameters['allLevelNames'] = [experimentParameters['columnVariableName']]+conditionNames
            experimentParameters['numConditionLevelValues'] = numConditionLevelValues
            experimentParameters['numericLevels'] = numericlevels
            experimentParameters[''] = tiledLevels
            master.switch_frame(columnLevelValuesPage,folderName)

        def backCommand():
            if experimentParameters['format'] == 'tube':
                master.switch_frame(TubeExperimentParameterPage,folderName)
            else:
                master.switch_frame(PlateExperimentParameterPage,folderName)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=numAllLevels+1,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: backCommand()).grid(row=numAllLevels+1,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=numAllLevels+1,column=2)

class columnLevelValuesPage(tk.Frame):
    def __init__(self, master,folderName):
        numColumnLevelValues = experimentParameters['numColumnLevelValues']
        
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        lt = tk.Label(mainWindow,text=experimentParameters['columnVariableName']+':').grid(row=0,column=0)
        col_wrap = 12
        for col in range(1,numColumnLevelValues+1):
            lt1 = tk.Label(mainWindow, text='Level Value '+str(col),width=10).grid(row=int((col-1)/col_wrap)*2,column=((col-1)%col_wrap)+1)
        entryList = []
        if 'sampleNameFile.xlsx' in os.listdir('misc'):
            df = pd.read_excel('misc/sampleNameFile.xlsx')
            if 'Time' in df.columns:
                times = list(map(float,list(pd.unique(df['Time']))))
                sortedTimes = list(map(str,sorted(times)))
                j=0

        for columnLevelValueNumber in range(numColumnLevelValues):
            e1 = tk.Entry(mainWindow,width=10)
            e1.grid(row=2*int(columnLevelValueNumber/col_wrap)+1,column=(columnLevelValueNumber%col_wrap)+1)
            if 'sampleNameFile.xlsx' in os.listdir('misc'):
                if 'Time' in df.columns:
                    e1.insert(tk.END,sortedTimes[j]) 
                    j+=1
            entryList.append(e1)

        def collectInputs():
            columnLevelValues = []
            for entry in entryList:
                columnLevelValues.append(float(entry.get()))
            experimentParameters['columnLevelValues'] = columnLevelValues
            master.switch_frame(conditionLevelValuesPage,folderName)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=2*int(numColumnLevelValues/col_wrap)+2,column=5)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(allLevelNamePage,folderName)).grid(row=2*int(numColumnLevelValues/col_wrap)+2,column=6)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=2*int(numColumnLevelValues/col_wrap)+2,column=7)

class LinkedEntryButton(tk.Button):
    def __init__(self,parent,linkedentries,**kwargs):
        tk.Button.__init__(self,parent,**kwargs)
        self.parent = parent
        self.linkedentries = linkedentries

    def fillNumerically(self):
        for i,entry in enumerate(self.linkedentries):
            entry.insert(tk.END,str(i+1))

class conditionLevelValuesPage(tk.Frame):
    def __init__(self, master,folderName):
        numConditionLevels = experimentParameters['numConditionLevels']
        maxLevelValues = max(experimentParameters['numConditionLevelValues'])
        
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        entryWrap = 16
        fullEntryList = []
        blist = []

        labelList = experimentParameters['conditionLevelNames']
        rowNum = 0
        if 'sampleNameFile.xlsx' in os.listdir('misc'):
            df = pd.read_excel('misc/sampleNameFile.xlsx')
        for conditionLevelNumber in range(numConditionLevels):
            if 'sampleNameFile.xlsx' in os.listdir('misc'):
                j = 0
            currentLabel = labelList[conditionLevelNumber]
            initialRowNum = rowNum
            l1 = tk.Label(mainWindow, text='Level values for \"'+currentLabel+'\":').grid(row=rowNum,column=1,sticky=tk.W) 
            levelEntryList = []
            rowNum-=1
            for col in range(1,maxLevelValues+1):
                if 'sampleNameFile.xlsx' in os.listdir('misc'):
                    for column in df.columns:
                        if column == currentLabel:
                            filledInLevelValues = list(pd.unique(df[column]))
                if (col-1)%entryWrap == 0:
                    rowNum+=1
                if col < experimentParameters['numConditionLevelValues'][conditionLevelNumber]+1:
                    e1 = tk.Entry(mainWindow,width=8)
                    e1.grid(row=rowNum,column=(col-1)%entryWrap+2,sticky=tk.W)
                    if 'sampleNameFile.xlsx' in os.listdir('misc'):
                        e1.insert(tk.END,filledInLevelValues[j])
                        j+=1
                    levelEntryList.append(e1)
            b = LinkedEntryButton(mainWindow,levelEntryList,text='Fill Values Numerically')
            b.configure(command=b.fillNumerically)
            b.grid(row=initialRowNum,column=0,sticky=tk.W)
            blist.append(b)
            fullEntryList.append(levelEntryList)
            rowNum+=1

        def collectInputs():
            conditionLevels = {}
            for lvlentrylist,i in zip(fullEntryList,range(numConditionLevels)):
                tempLevels = []
                for entry in lvlentrylist:
                    tempLevels.append(str(entry.get()))
                conditionLevels[experimentParameters['conditionLevelNames'][i]] = tempLevels
            experimentParameters['conditionLevelValues'] = conditionLevels.copy()
            experimentParameters['allLevelValues'] = conditionLevels.copy()
            experimentParameters['allLevelValues'][experimentParameters['columnVariableName']] = experimentParameters['columnLevelValues']
            shortenedExperimentParameters = {}
            for keyToKeep in ['format','numPlates','overallPlateDimensions','barcodingDict','unpackingDict']:
                if keyToKeep in list(experimentParameters.keys()):
                    shortenedExperimentParameters[keyToKeep] = experimentParameters[keyToKeep]
            shortenedExperimentParameters['levelLabelDict'] = experimentParameters['allLevelValues']
            for dataType in dataTypeList:
                with open('misc/experimentParameters-'+folderName+'-'+dataType+'.json', 'w') as fp:
                    json.dump(shortenedExperimentParameters, fp)
            #{"A1-12": {"A1-4": ["Barcode 1-", "Barcode 2-"], "A5-8": ["Barcode 1+", "Barcode 2-"], "A9-12": ["Barcode 1-", "Barcode 2+"]}}
            #{"A1-4": ["A1", "A2", "A4", "A3"], "A5-8": ["A5", "A6", "A8", "A7"], "A9-12": ["A9", "A10", "A12", "A11"]}
            if 'cell' in dataTypeList and experimentParameters['format'] == 'plate':
                #Both
                if '+' in globalMultiplexingVar:
                    for plate in experimentParameters['barcodingDict']:
                        if plate not in os.listdir('inputData/singleCellCSVFiles/'):
                            subprocess.run(['mkdir','inputData/singleCellCSVFiles/'+plate])
                else:
                    if 'Barcoding' in globalMultiplexingVar:
                        for plate in experimentParameters['barcodingDict']:
                            if plate not in os.listdir('inputData/singleCellCSVFiles/'):
                                subprocess.run(['mkdir','inputData/singleCellCSVFiles/'+plate])
                    elif '96' in globalMultiplexingVar:
                        for plate in experimentParameters['unpackingDict']:
                            if plate not in os.listdir('inputData/singleCellCSVFiles/'):
                                subprocess.run(['mkdir','inputData/singleCellCSVFiles/'+plate])
                    else:
                        for plateNum in range(1,experimentParameters['numPlates']+1):
                            if 'A'+str(plateNum) not in os.listdir('inputData/singleCellCSVFiles/'):
                                subprocess.run(['mkdir','inputData/singleCellCSVFiles/A'+str(plateNum)])
            master.switch_frame(ExperimentSetupStartPage,folderName,backPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=numConditionLevels+1,column=int(maxLevelValues/2))
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(columnLevelValuesPage,folderName)).grid(row=numConditionLevels+1,column=int(maxLevelValues/2)+1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=numConditionLevels+1,column=int(maxLevelValues/2)+2)
