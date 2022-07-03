#!/usr/bin/env python3  
import os,pickle,shutil
import statistics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#from logicle import logicle,quantile 
from datetime import datetime
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import matplotlib.lines as linesmpl
from plateypus.dataprocessing.miscFunctions import find_nearest,returnTicks,returnGates,returnGatesLinear
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
import os
if os.name == 'nt':
    dirSep = '\\'
else:
    dirSep = '/'

class ProliferationSelectionPage(tk.Frame):
    def __init__(self,master,fName,expn,exdata,shp,bp):
        tk.Frame.__init__(self, master)
        
        global secondaryhomepage
        secondaryhomepage = shp
        global folderName
        folderName = fName
        global expNum
        expNum = expn
        global ex_data
        ex_data = exdata
        global bPage
        bPage = bp
        
        singleCellDf = pd.read_hdf('outputData'+dirSep+'pickleFiles'+dirSep+'initialSingleCellDf-channel-'+folderName+'.h5',key='df')
        markers = list(singleCellDf.columns)
        celltypes = singleCellDf.index.unique('CellType').tolist()
        
        mainWindow = tk.Frame(self)
        mainWindow.pack()

        tk.Label(mainWindow,text='Proliferating:').grid(row=0,column=0,columnspan=2,sticky=tk.W)
        tk.Label(mainWindow,text='Celltypes:').grid(row=1,column=0,sticky=tk.W)
        rbList,rbVarList = [],[]
        for i,celltype in enumerate(celltypes):
            rbVar = tk.BooleanVar(value=False)
            rb = tk.Checkbutton(mainWindow,text=celltype,variable=rbVar)
            rb.grid(row=i+2,column=0,sticky=tk.W)
            rbList.append(rb)
            rbVarList.append(rbVar)
        
        tk.Label(mainWindow,text='Markers:').grid(row=1,column=1,sticky=tk.W)
        rbList2,rbVarList2 = [],[]
        for i,marker in enumerate(markers):
            rbVar = tk.BooleanVar(value=False)
            rb = tk.Checkbutton(mainWindow,text=marker,variable=rbVar)
            rb.grid(row=i+2,column=1,sticky=tk.W)
            rbList2.append(rb)
            rbVarList2.append(rbVar)

        def collectInputs():
            prolifCellTypes = [celltypes[i] for i,x in enumerate(rbVarList) if x.get()]
            prolifMarkers = [markers[i] for i,x in enumerate(rbVarList2) if x.get()]
            #Only select first prolif marker and first prolif celltype for now. TODO: allow multiple proliferation markers/celltypes
            global logicleDataStacked 
            logicleDataStacked = singleCellDf[prolifMarkers[0]].xs((prolifCellTypes[0]),level=('CellType'))
            master.switch_frame(GatingPage)
         

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(secondaryhomepage,folderName,expNum,ex_data,bPage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class GatingPage(tk.Frame):
    def __init__(self,master):
        tk.Frame.__init__(self, master)
        
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        v1 = tk.StringVar(value='Condition-Time')
        l1 = tk.Label(mainWindow,text = 'Do you want to group all samples, all timepoints, or all conditions together for initial CTV Gating?').grid(row=0,column=0,sticky=tk.W)
        rb1 = tk.Radiobutton(mainWindow,text='Group all Samples',variable=v1,value='Condition-Time')
        rb2 = tk.Radiobutton(mainWindow,text='Group all conditions',variable=v1,value='Condition')
        rb3 = tk.Radiobutton(mainWindow,text='Group all timepoints',variable=v1,value='Time')
        rb1.grid(row=1,column=0,sticky=tk.W)
        rb2.grid(row=2,column=0,sticky=tk.W)
        rb3.grid(row=3,column=0,sticky=tk.W)
        
        tk.Label(mainWindow,text='Maximum number of generations').grid(row=4,column=0)
        e1 = tk.Entry(mainWindow)
        e1.insert(tk.END,'6')
        e1.grid(row=5,column=0)

        #logicleDataStacked = pickle.load(open('outputData/pickleFiles/logicleProliferationDf.pkl','rb'))
        groupedGateList = []

        rowAsAllLevelNames = logicleDataStacked.index.names[:-1]
        rowAsTimePointLevelNames = logicleDataStacked.index.names[-2]
        rowAsConditionLevelNames = logicleDataStacked.index.names[:-2] 
        conditionTuple = []
        for tpl in logicleDataStacked.groupby(level=rowAsConditionLevelNames,sort=False):
            conditionTuple.append(tpl[0])
        timepointTuple = []
        for tpl in logicleDataStacked.groupby(level=rowAsTimePointLevelNames,sort=False):
            timepointTuple.append(tpl[0])
        
        #os.chdir('../../experiments/'+folderName) 
        def collectInputs():
            global groupVariable
            global titleVariable
            global outerVariableValues
            global innerVariableValues
            global titleVariable
            global levelNames
            global numGens
            groupVariable = v1.get()
            numGens = int(e1.get())
            #Group all samples together (no iteration at all)
            if groupVariable == 'Condition-Time':
                levelNames = rowAsAllLevelNames
                outerVariableValues = conditionTuple
                innerVariableValues = timepointTuple
                titleVariable = 'All'
                iterationRange = 1
            else:
                #Group conditions together (iterate through each timepoint):
                if groupVariable == 'Condition':
                    levelNames = rowAsTimePointLevelNames
                    outerVariableValues = timepointTuple
                    innerVariableValues = conditionTuple
                    titleVariable = 'Time'
                #Group timepoints together (iterate through each condition):
                else:
                    levelNames = rowAsConditionLevelNames
                    outerVariableValues = conditionTuple
                    innerVariableValues = timepointTuple
                    titleVariable = 'Condition'
                global logicleDataUnstacked
                logicleDataUnstacked = logicleDataStacked.groupby(level=levelNames,sort=False).first().to_frame('MFI')
                iterationRange = logicleDataUnstacked.shape[0] 
            with open('misc'+dirSep+'iterationRange.pkl','wb') as f:
                pickle.dump(iterationRange,f)
            with open('misc'+dirSep+'rowValIterationRange.pkl','wb') as f:
                pickle.dump(0,f)
            groupedGateList = []
            with open('misc'+dirSep+'groupedGateList.pkl','wb') as f:
                pickle.dump(groupedGateList,f)
            master.switch_frame(Unsplit_Proliferation_Gating_GUI)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ProliferationSelectionPage,folderName,expNum,ex_data,secondaryhomepage,bPage)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class Unsplit_Proliferation_Gating_GUI(tk.Frame):
    def __init__(self,master):
        self.root = master.root
        tk.Frame.__init__(self, master)
        
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()
        
        plotFrame = tk.Frame(self)
        plotFrame.pack()
        row = pickle.load(open('misc'+dirSep+'rowValIterationRange.pkl','rb'))
        groupedGateList = pickle.load(open('misc'+dirSep+'groupedGateList.pkl','rb'))
        if groupVariable == 'Condition-Time':
            currentDf = logicleDataStacked.to_frame('MFI')
        else:
            currentLevelValues = logicleDataUnstacked.iloc[row,:].name
            currentDf = logicleDataStacked.xs(currentLevelValues,level=levelNames).to_frame('MFI')
        plottingDf = currentDf.reset_index()
        if groupVariable == 'Condition-Time' or groupVariable == 'Condition':
            nonEventLevelNames = currentDf.index.names[:-1]
            groupingColumn = []
            for conditionGroupbyTuple in currentDf.groupby(level=nonEventLevelNames,sort=False):
                conditionGroupbyTupleNew = [str(i) for i in conditionGroupbyTuple[0]]
                conditionName = '-'.join(conditionGroupbyTupleNew)
                for event in range(conditionGroupbyTuple[1].shape[0]):
                    groupingColumn.append(conditionName)
            plottingDf[groupVariable] = groupingColumn
        
        #Plot facetgrid of histograms of ctv values of all times for current condition
        kind='kde'
        if groupVariable == 'Condition-Time':
            g = sns.displot(data=plottingDf,x='MFI',color='purple',height=6,kind=kind)
        elif groupVariable == 'Time':
            g = sns.displot(data=plottingDf,x='MFI',hue='Time',height=6,palette=sns.color_palette("Purples", len(pd.unique(plottingDf[groupVariable]))),kind=kind)
        else:
            g = sns.displot(data=plottingDf,x='MFI',hue='Condition',height=6,palette=sns.color_palette("Purples", len(pd.unique(plottingDf['Condition']))),kind=kind)
        #g = sns.FacetGrid(plottingDf,legend_out=True,hue=groupVariable,height=6,palette=sns.color_palette("Purples", len(pd.unique(plottingDf[groupVariable]))))
        #g.map(sns.distplot,'MFI',kde=False,bins=256)
        #g.map(sns.kdeplot,'MFI',bw=15)
        plt.subplots_adjust(top=0.95)
        xtickValues,xtickLabels = returnTicks([-1000,100,1000,10000,100000])
        axis = g.fig.get_axes()[0]
        axis.set_xticks(xtickValues)
        axis.set_xticklabels(xtickLabels)
        
        #startingGateLogicle = 725
        space = np.linspace(0,1024,num=100)
        kde = gaussian_kde(currentDf['MFI'].values)
        kdeVals = kde.evaluate(space)
        maxInd = np.argmax(kdeVals)
        minMaxScaledKDE = kdeVals.reshape(-1,1)
        minMaxScaledKDE = (minMaxScaledKDE - minMaxScaledKDE.min()) / (minMaxScaledKDE.max() - minMaxScaledKDE.min())
        #tempDf = pd.DataFrame({sizeDimension:space,'Density':minMaxScaledKDE[:,0]})
        startingGateLogicle = space[np.argmax(minMaxScaledKDE[:,0] > 0.01)]
        endingGateLogicle = max(space)

        logicleGates = [startingGateLogicle]+[space[x] for x in list(find_peaks(kdeVals*-1)[0]) if space[x] > startingGateLogicle]+[endingGateLogicle]
        logicleGates = logicleGates[::-1]
        #logicleGates = returnGatesLinear(currentDf,startingGateLogicle,numGens)
        maxY = list(plt.gca().get_yticks())[-1]
        X = endingGateLogicle
        Ymin = -1
        Ymax = maxY
        gates = logicleGates[1:]
        
        #sys.exit(0)
        numGates = len(logicleGates)-1
        self.canvas = FigureCanvasTkAgg(g.fig,master=plotFrame)
        self.ax = plt.gca() 
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()
        self.X = X
        self.gates= gates
        self.gateVals = gates
        gateLengths = [y-x for x, y in zip(gates[:-1], gates[1:])][1:]
        self.gateLengths = gateLengths 
        self.numGates = numGates
        x = [X, X]
        y = [Ymin, Ymax]
        
        self.line = linesmpl.Line2D(x, y, picker=5,linestyle=':',color='r')
        self.ax.add_line(self.line)
        self.childlines = []
        for gate in range(self.numGates):
            newX = self.gates[gate]
            childline = linesmpl.Line2D([newX,newX], y,linestyle=':',color='k',picker=5)
            self.childlines.append(childline)
        for childline in self.childlines:
            self.ax.add_line(childline)
        self.canvas.draw()
        self.sid = self.canvas.mpl_connect('pick_event', self.clickOnParentLine)
        childsids = []
        for childline,i in zip(self.childlines,range(self.numGates)):
            childsids.append(self.canvas.mpl_connect('pick_event', self.clickOnChildLines))
        self.childsids = childsids
        
        #Add title
        if groupVariable == 'Condition-Time':
            plt.title(titleVariable)
        else:
            if not isinstance(currentLevelValues, (tuple,)):
                currentLevelValues = [str(currentLevelValues)]
            plt.title('-'.join(currentLevelValues)+' ('+titleVariable+' '+str(row+1)+'/'+str(logicleDataUnstacked.shape[0])+')')
        
        def collectInputs():
            for names in range(len(pd.unique(plottingDf[groupVariable]))):
                groupedGateList.append(self.getAllX())
            with open('misc'+dirSep+'groupedGateList.pkl','wb') as f:
                pickle.dump(groupedGateList,f)
            row = pickle.load(open('misc'+dirSep+'rowValIterationRange.pkl','rb'))
            row+=1
            with open('misc'+dirSep+'rowValIterationRange.pkl','wb') as f:
                pickle.dump(row,f)
            if row != pickle.load(open('misc'+dirSep+'iterationRange.pkl','rb')):
                master.switch_frame(Unsplit_Proliferation_Gating_GUI)
            else:
                i = 0
                row = 0
                singleCellProliferationDf = pd.DataFrame(np.zeros(logicleDataStacked.shape),logicleDataStacked.index,columns=['Generation'])
                for outerVariable in outerVariableValues:
                    for innerVariable in innerVariableValues:
                        sampleGenerationGates = groupedGateList[i]
                        #sampleGenerationGates = groupedGate.getAllX()
                        if groupVariable == 'Condition':
                            indexingTuple = tuple(list(innerVariable)+[outerVariable,slice(None)])
                        else:
                            indexingTuple = tuple(list(outerVariable)+[innerVariable,slice(None)])
                        ctvValues = logicleDataStacked.loc[indexingTuple]
                        generationValues = np.zeros(ctvValues.shape)
                        for sampleEvent,row in zip(ctvValues,range(ctvValues.shape[0])):
                            for generation in range(len(sampleGenerationGates)-1):
                                upperGate = sampleGenerationGates[generation]
                                lowerGate = sampleGenerationGates[generation+1]
                                if(sampleEvent > lowerGate and sampleEvent <= upperGate):
                                    generationValues[row] = generation
                        singleCellProliferationDf.loc[indexingTuple,:] = generationValues

                print(singleCellProliferationDf)
                with open('outputData'+dirSep+'pickleFiles'+dirSep+'singleCellDataFrame-proliferation-'+folderName+'.pkl', 'wb') as f:
                    pickle.dump(singleCellProliferationDf,f)
                master.switch_frame(secondaryhomepage,folderName,expNum,ex_data,bPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

    def clickOnParentLine(self, event):
        if event.artist == self.line:
            #print("line selected ", event.artist)
            self.follower = self.canvas.mpl_connect("motion_notify_event", self.followParentMouse)
            self.releaser = self.canvas.mpl_connect("button_press_event", self.releaseParentOnClick)
    
    def clickOnChildLines(self, event):
        if event.artist in self.childlines:
            #print("line selected ", event.artist)
            self.currentArtist = event.artist
            self.follower = self.canvas.mpl_connect("motion_notify_event", self.followChildMouse)
            self.releaser = self.canvas.mpl_connect("button_press_event", self.releaseChildOnClick)
    
    def followChildMouse(self, event):
        for childline,gate in zip(self.childlines,range(self.numGates)):
            if self.childlines[gate] == self.currentArtist:
                self.childlines[gate].set_xdata([event.xdata, event.xdata])
        self.canvas.draw()
        #self.c.draw_idle()
    
    def releaseChildOnClick(self, event):
        newGates = []
        for childline in self.childlines:
            newGates.append(childline.get_xdata()[0])
        self.gates = newGates
        self.canvas.mpl_disconnect(self.releaser)
        self.canvas.mpl_disconnect(self.follower)

    def followParentMouse(self, event):
        self.line.set_xdata([event.xdata, event.xdata])
        newGates = []
        for childline,gate in zip(self.childlines,range(self.numGates)):
            newX = (event.xdata-self.X)+self.gates[gate]
            childline.set_xdata([newX,newX])
            newGates.append(newX)
        self.gateVals = newGates
        self.canvas.draw()
    
    def releaseParentOnClick(self, event):
        self.X = self.line.get_xdata()[0]
        newGates = []
        for childline in self.childlines:
            newGates.append(childline.get_xdata()[0])
        self.gates = newGates
        self.canvas.mpl_disconnect(self.releaser)
        self.canvas.mpl_disconnect(self.follower)

    def getParentX(self):
        return self.X
    
    def getChildX(self):
        childGates = []
        for childline in self.childlines:
            childGates.append(childline.get_xdata()[0])
        return childGates

    def getAllX(self):
        return [self.getParentX()]+self.getChildX()

def generateBulkProliferationStatistics(folderName,experimentNumber):
    singleCellProliferationDfStacked = pickle.load(open('outputData'+dirSep+'pickleFiles'+dirSep+'singleCellDataFrame-proliferation-'+folderName+'.pkl','rb'))
    idx = pd.IndexSlice
    numGenerations = len(singleCellProliferationDfStacked.loc[:,'Generation'].unique())
    singleCellProliferationDfUnstacked = singleCellProliferationDfStacked.unstack('Event')
    #Reindex unstacked dataframe by original ordering
    #conditionLevelNames = singleCellProliferationDfStacked.index.names[:-1]
    #indexDf = singleCellProliferationDfStacked.groupby(level=conditionLevelNames,sort=False).first()
    #print(indexDf)
    #indexDfFull = pickle.load(open('outputData/pickleFiles/cellStatisticPickleFile-'+folderName+'.pkl','rb'))
    firstTimepoint = singleCellProliferationDfStacked.index.unique('Time').tolist()[0]
    reindex = singleCellProliferationDfStacked.xs((1,firstTimepoint),level=('Event','Time')).index
    indexDf = singleCellProliferationDfStacked.xs((1),level=('Event')).unstack('Time').reindex(reindex).loc[:,'Generation']
    #indexDf = indexDfFull.loc[list(pd.unique(indexDfFull.index.get_level_values(0)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(1)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(2)))[0]]
    levelValues = []
    iteration = 0
    #temp = indexDfFull.reset_index()
    temp = singleCellProliferationDfUnstacked.reset_index()
    emptyMatrixLogicle = np.zeros(singleCellProliferationDfUnstacked.shape)
    for row in range(indexDf.shape[0]):
        levelName = list(indexDf.iloc[row,:].name)
        for timepoint in list(pd.unique(indexDf.columns)):
            tp = tuple(levelName+[timepoint,slice(None)])
            emptyMatrixLogicle[iteration,:] = singleCellProliferationDfUnstacked.loc[tp,slice(None)]
            levelValues.append(tp[:-1])
            iteration+=1
    multiIndex = pd.MultiIndex.from_tuples(levelValues, names=indexDf.index.names+[indexDf.columns.name])
    
    singleCellProliferationDfUnstacked = pd.DataFrame(emptyMatrixLogicle,index=multiIndex,columns=singleCellProliferationDfUnstacked.columns)
    singleCellProliferationDfUnstacked.columns.name = 'Event'
    #Frequency
    generationFrequencyMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],numGenerations])
    for row in range(singleCellProliferationDfUnstacked.shape[0]):
        proliferationDfPerCondition = singleCellProliferationDfUnstacked.iloc[row,:]
        freqTable = pd.value_counts(proliferationDfPerCondition.values.ravel()).sort_index()
        generations = [float(i) for i in list(freqTable.index.get_level_values(0))]
        for generation in generations:
            generationFrequencyMatrix[row,int(generation)] = freqTable[generation]
    generationFrequencyDataFrame = pd.DataFrame(generationFrequencyMatrix,index=singleCellProliferationDfUnstacked.index,columns=range(numGenerations)).astype(int)
    generationFrequencyDataFrame.columns.name = 'Generation'
    
    fractionDilutedMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],1])
    precursorFrequencyMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],1])
    proliferationIndexMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],1])
    divisionIndexMatrix = np.zeros([singleCellProliferationDfUnstacked.shape[0],1])
    for row in range(generationFrequencyDataFrame.shape[0]):
        precursorFrequencyNumerator = 0
        proliferationIndexNumerator = 0
        for generation in range(1,numGenerations):
            precursorFrequencyNumerator+=generationFrequencyDataFrame.iloc[row,generation]/(2**generation)
            proliferationIndexNumerator+=generation*generationFrequencyDataFrame.iloc[row,generation]/(2**generation)
        precursorFrequencyDenominator = 0
        for generation in range(numGenerations):
            precursorFrequencyDenominator+=generationFrequencyDataFrame.iloc[row,generation]/(2**generation)
        
        fractionDiluted = np.sum(generationFrequencyDataFrame.iloc[row,1:])/np.sum(generationFrequencyDataFrame.iloc[row,0:])
        precursorFrequency = precursorFrequencyNumerator/precursorFrequencyDenominator
        proliferationIndex = proliferationIndexNumerator/precursorFrequencyNumerator
        divisionIndex = proliferationIndexNumerator/precursorFrequencyDenominator

        fractionDilutedMatrix[row,0] = fractionDiluted
        precursorFrequencyMatrix[row,0] = precursorFrequency
        proliferationIndexMatrix[row,0] = proliferationIndex
        divisionIndexMatrix[row,0] = divisionIndex
    allProliferationStatisticsList = []
    allProliferationMatricesList = [fractionDilutedMatrix,precursorFrequencyMatrix,proliferationIndexMatrix,divisionIndexMatrix]
    allProliferationStatisticNames = ['fractionDiluted','precursorFrequency','proliferationIndex','divisionIndex']
    for i in range(generationFrequencyDataFrame.shape[1]):
        currentGenerationFrequency = np.array(generationFrequencyDataFrame.iloc[:,i].values.ravel())
        allProliferationStatisticNames.append('rawFrequencyGen'+str(i))
        allProliferationMatricesList.append(currentGenerationFrequency)
        allProliferationStatisticNames.append('normalizedFrequencyGen'+str(i))
        allProliferationMatricesList.append(np.divide(currentGenerationFrequency,2**i))
    for i in range(len(allProliferationMatricesList)):
        statisticDf = pd.DataFrame(allProliferationMatricesList[i],index=singleCellProliferationDfUnstacked.index)#,columns = [allProliferationStatisticNames[i]])
        statisticDfUnstacked = statisticDf.unstack('Time')
        #indexDfFull = pickle.load(open('outputData/pickleFiles/cellStatisticPickleFile-'+folderName+'.pkl','rb'))
        #indexDf = indexDfFull.loc[list(pd.unique(indexDfFull.index.get_level_values(0)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(1)))[0]].loc[list(pd.unique(indexDfFull.index.get_level_values(2)))[0]]
        levelValues = []
        iteration = 0
        emptyMatrixLogicle = np.zeros(statisticDfUnstacked.shape)
        for row in range(indexDf.shape[0]):
            levelName = list(indexDf.iloc[row,:].name)
            tp = tuple(levelName)
            emptyMatrixLogicle[iteration,:] = statisticDfUnstacked.loc[tp,:]
            levelValues.append(tp)
            iteration+=1
        multiIndex = pd.MultiIndex.from_tuples(levelValues, names=indexDf.index.names)
        
        statisticDfUnstacked = pd.DataFrame(emptyMatrixLogicle,index=multiIndex,columns=statisticDfUnstacked.columns)
        statisticDfUnstacked.columns.name = 'Time'
        allProliferationStatisticsList.append(statisticDfUnstacked)
    completeBulkProliferationDataFrame = pd.concat(allProliferationStatisticsList,keys=allProliferationStatisticNames,names=['Statistic'],axis=0)
    completeBulkProliferationDataFrame.columns = completeBulkProliferationDataFrame.columns.droplevel(0)
    return completeBulkProliferationDataFrame
