#!/usr/bin/env python3
import os,pickle
import pandas as pd
import numpy as np
import tkinter as tk
from plateypus.dataprocessing.miscFunctions import setMaxWidth

expParamDict = {'cyt':'cyt','cell':'cell','prolif':'cell','singlecell':'cell'}

class TubeLayoutPage(tk.Frame):
    def __init__(self,master,folderName,levelLabelDict,numSamples,dataType,secondaryhomepage,backPage):
        allLevelNames = list(levelLabelDict.keys())
        conditionLevelValues = [levelLabelDict[x] for x in list(levelLabelDict.keys())[:-1]]
        columnLevelValues = levelLabelDict[list(levelLabelDict.keys())[-1]]
        numericList = [False]*len(conditionLevelValues) + [True]

        tk.Frame.__init__(self, master)
        
        labelWindow = tk.Frame(self)
        labelWindow.pack()
        l1 = tk.Label(labelWindow, text="Assign level values to each sample:", font='Helvetica 18 bold').pack()

        """BEGIN TEMP SCROLLBAR CODE"""
        labelWindow1 = tk.Frame(self)
        labelWindow1.pack(side=tk.TOP,padx=10,fill=tk.X,expand=True)
        
        #Make canvas
        w1 = tk.Canvas(labelWindow1, width=1500, height=400,background="white", scrollregion=(0,0,3000,33*numSamples))

        #Make scrollbar
        scr_v1 = tk.Scrollbar(labelWindow1,orient=tk.VERTICAL)
        scr_v1.pack(side=tk.RIGHT,fill=tk.Y)
        scr_v1.config(command=w1.yview)
        #Add scrollbar to canvas
        w1.config(yscrollcommand=scr_v1.set)
        
        scr_v2 = tk.Scrollbar(labelWindow1,orient=tk.HORIZONTAL)
        scr_v2.pack(side=tk.BOTTOM,fill=tk.X)
        scr_v2.config(command=w1.xview)
        w1.config(xscrollcommand=scr_v2.set)

        w1.pack(fill=tk.BOTH,expand=True)
        #Make and add frame for widgets inside of canvas
        #canvas_frame = tk.Frame(w1)
        mainWindow = tk.Frame(w1)
        mainWindow.pack()
        w1.create_window((0,0),window=mainWindow, anchor = tk.NW)
        """END TEMP SCROLLBAR CODE"""
 
        #mainWindow = tk.Frame(self)
        #mainWindow.pack(fill=tk.BOTH, expand=True)
        #scrollable_body = Scrollable(mainWindow, width=100)
        
        #Automatically import file names into layout dict if excel/csv template in correct folder
        if 'sampleNameFile.xlsx' in os.listdir('misc') or 'sampleNameFile.csv' in os.listdir('misc'):
            if 'sampleNameFile.xlsx' in os.listdir('misc'): 
                sampleNameDf = pd.read_excel('misc/sampleNameFile.xlsx')
            else:
                sampleNameDf = pd.read_csv('misc/sampleNameFile.csv')
        else:
            sampleNameDf = []

        #Should be in this format: headers are: sampleFileName Level0Name Level1Name
        fcsFiles = []
        if 'A1_cell.csv' not in os.listdir('inputData/bulkCSVFiles/'):
            for fcsName in os.listdir('inputData/fcsFiles'):
                if '.fcs' in fcsName:
                    fcsFiles.append(fcsName)
            if len(fcsFiles) == 0:
                fcsFiles = ['                 ']
        else:
            bulkStatFile = pd.read_csv('inputData/bulkCSVFiles/A1_cell.csv')
            for row in range(bulkStatFile.shape[0]):
                if bulkStatFile.iloc[row,0] not in ['Mean','SD']:
                    fcsFiles.append(bulkStatFile.iloc[row,0])

        fcsMenuList = []
        fcsVarList = []
        includeCheckboxList = []
        includeBoolVarList = []
        sampleDropdownList = []
        sampleParameterList = []
        
        tk.Label(mainWindow,text='FileName: ').grid(row=0,column=1)
        tk.Label(mainWindow,text='Include?').grid(row=0,column=2)
        for i,level in enumerate(allLevelNames):
            tk.Label(mainWindow,text=level+': ').grid(row=0,column=i+3)
            
        for sample in range(numSamples):
            tk.Label(mainWindow,text='Sample '+str(sample+1)+': ').grid(row=sample+1,column=0)
            fcsVar = tk.StringVar()
            #fcsMenu = tk.OptionMenu(scrollable_body,fcsVar,*fcsFiles)
            fcsMenu = tk.OptionMenu(mainWindow,fcsVar,*fcsFiles)
            fcsMenu.grid(row=sample+1,column=1)
            setMaxWidth(fcsFiles,fcsMenu)
            
            fcsMenuList.append(fcsMenu)
            fcsVarList.append(fcsVar)
            includeBoolVar = tk.BooleanVar(value=True)
            includeCheckbutton = tk.Checkbutton(mainWindow,variable=includeBoolVar)
            includeCheckbutton.grid(row=sample+1,column=2)
            includeCheckboxList.append(includeCheckbutton)
            includeBoolVarList.append(includeBoolVar)
            
            currentSampleDropdownList = []
            currentSampleParameterList = []
            for i,level in enumerate(allLevelNames):
                if level != 'Time':
                    levelValues = conditionLevelValues[i]
                else:
                    levelValues = columnLevelValues
                levelValueVar = tk.StringVar()
                #levelValueMenu = tk.OptionMenu(scrollable_body,levelValueVar,*levelValues)
                levelValueMenu = tk.OptionMenu(mainWindow,levelValueVar,*levelValues)
                if len(levelValues) == 1:
                    levelValueVar.set(levelValues[0])
                else:
                    if not isinstance(sampleNameDf,list) and level in sampleNameDf.columns:
                       fcsVar.set(sampleNameDf.iloc[sample,0])
                       levelValueVar.set(sampleNameDf.iloc[sample,:][level])
                levelValueMenu.grid(row=sample+1,column=i+3)
                setMaxWidth(levelValues,levelValueMenu)
                currentSampleDropdownList.append(levelValueMenu)
                currentSampleParameterList.append(levelValueVar)
            sampleDropdownList.append(currentSampleDropdownList)
            sampleParameterList.append(currentSampleParameterList)

        #scrollable_body.update() 
        
        def collectInputs():
            fullSampleList = []
            indexTuples = []
            for i,sample in enumerate(fcsVarList):
                if includeBoolVarList[i].get():
                    fullSampleList.append(sample.get())
                    sampleTuple = []
                    for j,sampleParameter in enumerate(sampleParameterList[i]):
                        if numericList[j]:
                            sampleTuple.append(float(sampleParameter.get()))
                        else:
                            sampleTuple.append(sampleParameter.get())
                    indexTuples.append(sampleTuple)
            mi = pd.MultiIndex.from_tuples(indexTuples,names=allLevelNames)
            data = fullSampleList 
            df1 = pd.Series(data,index=mi)
            df1.to_pickle('misc/tempDf.pkl')
            df2 = df1.unstack('Time')
            
            #Reindex df after unstacking by original order
            xsdf = df1.xs([columnLevelValues[0]],level=['Time'])
            emptydf = np.empty(df2.shape,dtype='object')
            for row in range(emptydf.shape[0]):
                for col in range(emptydf.shape[1]):
                    timepoint = df2.columns[col]
                    if not isinstance(xsdf.index.tolist()[row],str) and not isinstance(xsdf.index.tolist()[row],float):
                        conditions = list(xsdf.index.tolist()[row])
                    else:
                        conditions = [xsdf.index.tolist()[row]]
                    emptydf[row,col] = df2.loc[tuple(conditions),timepoint]
            reindexedDf = pd.DataFrame(emptydf,index=xsdf.index,columns=df2.columns)
            reindexedDf.columns.name = 'Time'
            
            if dataType == 'both':
                dataTypes = ['cell','cyt']
            else:
                dataTypes = [dataType]
            for dt in dataTypes:
                with open('misc/tubeLayout-'+folderName+'-'+dt+'.pkl','wb') as f:
                    pickle.dump(reindexedDf,f)
            
            master.switch_frame(secondaryhomepage,folderName,backPage)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(secondaryhomepage,folderName,backPage)).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)
