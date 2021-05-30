#!/usr/bin/env python3
import pickle,os,json
import pandas as pd
import tkinter as tk
import tkinter.ttk
import plateypus.dataprocessing.initialDataProcessing as idp
import plateypus.dataprocessing.miscFunctions as mf
import plateypus.plotting.facetPlotLibrary as fpl 
import plateypus.plotting.interactiveGUIElements as ipe

expParamDict = {'cyt':'cyt','cell':'cell','prolif':'cell','singlecell':'cell','killing':'killing'}

#Get level names and values into an easily accessible dictionary
def createLabelDict(df):
    fulldf = df.stack()
    labelDict = {}
    for i in range(fulldf.index.nlevels):
        levelName = fulldf.index.levels[i].name
        if levelName not in ['Event','event']:
            labelDict[levelName] = list(pd.unique(fulldf.index.get_level_values(levelName)))
    return labelDict

def createLabelDictWithExperimentParameters(df,experimentParameters):
    fulldf = df.stack()
    labelDict = {}
    for i in range(fulldf.index.nlevels):
        levelName = fulldf.index.levels[i].name
        if levelName in ['Event','event']:
            pass
        else:
            if 'allLevelValues' in list(experimentParameters.keys()):
                experimentParameters['levelLabelDict'] = experimentParameters['allLevelValues']
            if levelName in experimentParameters['levelLabelDict'].keys():
                if len(experimentParameters['levelLabelDict'][levelName]) > 0:
                    labelDict[levelName] = experimentParameters['levelLabelDict'][levelName]
                else:
                    labelDict[levelName] = list(pd.unique(fulldf.index.get_level_values(levelName)))
            else:
                labelDict[levelName] = list(pd.unique(fulldf.index.get_level_values(levelName)))
    return labelDict

class checkUncheckAllButton(tk.Button):
    def __init__(self,parent,checkButtonList,**kwargs):
        tk.Button.__init__(self,parent,**kwargs)
        self.checkButtonList = checkButtonList
        self.parent = parent

    def checkAll(self):
        for checkButton in self.checkButtonList:
            checkButton.select()
    
    def uncheckAll(self):
        for checkButton in self.checkButtonList:
            checkButton.deselect()

class PlotExperimentWindow(tk.Frame):
    def __init__(self, master,fName,sPage):
        with open('misc/normalPlottingBool.pkl','wb') as f:
            pickle.dump(True,f)

        global folderName,switchPage
            
        folderName = fName
        switchPage = sPage

        tk.Frame.__init__(self, master)
        
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l2 = tk.Label(mainWindow, text="""Datatype: """)
        v2 = tk.StringVar(value='cyt')
        rb2a = tk.Radiobutton(mainWindow, text="Cytokine",padx = 20, variable=v2, value='cyt')
        rb2b = tk.Radiobutton(mainWindow,text="Cell",padx = 20, variable=v2, value='cell')
        rb2c = tk.Radiobutton(mainWindow,text="Proliferation",padx = 20, variable=v2, value='prolif')
        rb2d = tk.Radiobutton(mainWindow,text="Single Cell",padx = 20, variable=v2, value='singlecell')
        rb2e = tk.Radiobutton(mainWindow,text="Killing",padx = 20, variable=v2, value='killing')
        
        l2.grid(row=0,column=0)
        rb2a.grid(row=1,column=0,sticky=tk.W)
        rb2b.grid(row=2,column=0,sticky=tk.W)
        rb2c.grid(row=3,column=0,sticky=tk.W)
        rb2d.grid(row=4,column=0,sticky=tk.W)
        rb2e.grid(row=5,column=0,sticky=tk.W)
        
        def collectInputs():
            global useModifiedDf
            modified = False 
            useModifiedDf = modified
            if modified:
                modifiedString = '-modified'
            else:
                modifiedString = ''
            global dataType
            dataType = str(v2.get())
            global experimentDf
            if dataType != 'singlecell':
                experimentDf = pd.read_pickle('outputData/pickleFiles/'+idp.dataTypeDataFrameFileNames[dataType]+'-'+folderName+modifiedString+'.pkl')
            else:
                if 'initialSingleCellDf-channel-'+folderName+modifiedString+'.pkl' in os.listdir('outputData/pickleFiles/'):
                    experimentDf = pickle.load(open('outputData/pickleFiles/'+'initialSingleCellDf-channel-'+folderName+modifiedString+'.pkl','rb'))
                else:
                    experimentDf = pd.read_hdf('outputData/pickleFiles/'+'initialSingleCellDf-channel-'+folderName+modifiedString+'.h5', 'df')
                if 'CellType' != experimentDf.index.names[0]:
                    experimentDf = pd.concat([experimentDf], keys=['TCells'], names=['CellType'])
            global trueLabelDict
            trueLabelDict = {}
            if experimentDf.columns.name == None:
                experimentDf.columns.name = 'Marker'
            experimentParametersBool = False
            for fn in os.listdir('misc'):
                if 'experimentParameters' in fn:
                    if expParamDict[dataType] in fn:
                        experimentParametersBool = True
                        experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'-'+expParamDict[dataType]+'.json','r'))
            if experimentParametersBool:
                trueLabelDict = createLabelDictWithExperimentParameters(experimentDf,experimentParameters)
            else:
                trueLabelDict = createLabelDict(experimentDf)
            master.switch_frame(PlotTypePage)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(switchPage,folderName)).grid(row=5,column=1)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)


class PlotTypePage(tk.Frame):
    def __init__(self, master):

        if 'normalPlottingBool.pkl' in os.listdir('misc'):
            normalPlottingBool = pickle.load(open('misc/normalPlottingBool.pkl','rb'))
        else:
            normalPlottingBool = True

        if not normalPlottingBool:
            global useModifiedDf,dataType,experimentDf,trueLabelDict,folderName,switchPage
            plottingParams = pickle.load(open('misc/plottingParams.pkl','rb'))
            useModifiedDf = False
            dataType = 'cell'
            experimentDf = plottingParams['df']
            trueLabelDict = createLabelDict(experimentDf)
            folderName = plottingParams['folderName']
            switchPage = 'a' 

        plottableFigureDict = {'1d':['histogram','kde'],'categorical':['bar','violin','box','point','swarm','strip'],'2d':['line','scatter'],'3d':['heatmap']}
        
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        l1 = tk.Label(mainWindow, text="What type of figure do you want to plot?",pady=10).grid(row=0,column = 0,columnspan = len(plottableFigureDict),sticky=tk.N)
        
        #global trueLabelDict
        #trueLabelDict = {}
        #trueLabelDict = createLabelDict(experimentDf)
         
        plotTypeRadioButtons = []
        plotSelectionString = tk.StringVar(value='1d/histogram')
        maxNumSubplots = 0
        for plotTypeIndex,plotTypeTitle in enumerate(plottableFigureDict):
            plotTypeTitleLabel = tk.Label(mainWindow,text=plotTypeTitle).grid(row=1,column=plotTypeIndex,sticky=tk.NW)
            temprblist = []
            tempselectionstring = []
            for subPlotTypeIndex,subPlotTitle in enumerate(plottableFigureDict[plotTypeTitle]):
                rb = tk.Radiobutton(mainWindow, text=subPlotTitle,padx = 20, variable=plotSelectionString, value=plotTypeTitle+'/'+subPlotTitle)
                rb.grid(row=subPlotTypeIndex+2,column=plotTypeIndex,sticky=tk.NW)
                temprblist.append(rb)
            plotTypeRadioButtons.append(temprblist)
            if len(plottableFigureDict[plotTypeTitle]) > maxNumSubplots:
                maxNumSubplots = len(plottableFigureDict[plotTypeTitle])
        
        stripSwarmBool = tk.BooleanVar()
        cb = tk.Checkbutton(mainWindow,text='Add strip/swarm points to categorical plot',variable=stripSwarmBool,pady=20)
        cb.grid(row=maxNumSubplots+2,column=0,columnspan=len(plottableFigureDict))
        
        def collectInputs():
            global plotType
            global subPlotType
            global addDistributionPoints
            addDistributionPoints = stripSwarmBool.get()
            plotType,subPlotType = plotSelectionString.get().split('/')
            master.switch_frame(selectLevelsPage,'a','b','c','d','e')
        
        def backCommand():
            #os.chdir('../../')
            if normalPlottingBool:
                master.switch_frame(PlotExperimentWindow,folderName,switchPage)
            else:
                master.switch_frame(plottingParams['homepage'],folderName,plottingParams['bp'],plottingParams['shp'])
            #master.switch_frame(PlotExperimentWindow,switchPage)
       
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: backCommand()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)

class selectLevelsPage(tk.Frame):
    def __init__(self, master,fsp,fName,backPage,pt,shp):
        tk.Frame.__init__(self, master)
        labelWindow = tk.Frame(self)
        labelWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        global figureLevelList,fullFigureLevelBooleanList
        fullFigureLevelBooleanList = []
        figureLevelList = []
        
        l1 = tk.Label(labelWindow, text="""Which levels names do you want to be included within this figure??:""").pack()
        mainWindow = tk.Frame(self)
        levelNameCheckButtons = []
        checkButtonVariableList = []
        for levelName,i in zip(trueLabelDict.keys(),range(len(trueLabelDict.keys()))):
            includeLevelBool = tk.BooleanVar(value=True)
            cb = tk.Checkbutton(mainWindow, text=levelName,padx = 20, variable=includeLevelBool,onvalue=True)
            cb.grid(row=i+3,column=1,sticky=tk.W)
            cb.select()
            levelNameCheckButtons.append(cb)
            checkButtonVariableList.append(includeLevelBool)
        
        checkButtonWindow = tk.Frame(self)
        checkAllButton1 = checkUncheckAllButton(checkButtonWindow,levelNameCheckButtons, text='Check All')
        checkAllButton1.configure(command=checkAllButton1.checkAll)
        checkAllButton1.pack(side=tk.LEFT)
        
        uncheckAllButton1 = checkUncheckAllButton(checkButtonWindow,levelNameCheckButtons, text='Uncheck All')
        uncheckAllButton1.configure(command=checkAllButton1.uncheckAll)
        uncheckAllButton1.pack(side=tk.LEFT)
        
        checkButtonWindow.pack(side=tk.TOP)
        mainWindow.pack(side=tk.TOP,padx=10)

        def collectInputs():
            includeLevelList = []
            for checkButtonVariable in checkButtonVariableList:
                includeLevelList.append(checkButtonVariable.get())
            for figureLevelBool,levelName in zip(includeLevelList,trueLabelDict):
                if figureLevelBool:
                    figureLevelList.append(levelName)
                    fullFigureLevelBooleanList.append(True)
                else:
                    fullFigureLevelBooleanList.append(False)
            master.switch_frame(selectLevelValuesPage,assignLevelsToParametersPage,trueLabelDict,selectLevelsPage,selectLevelsPage,fsp,fName,shp,'a')
        
        def quitCommand():
            exitBoolean = True
            quit()

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(PlotTypePage)).pack(in_=buttonWindow,side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quitCommand()).pack(in_=buttonWindow,side=tk.LEFT)

class selectLevelValuesPage(tk.Frame):
    #(selectLevelValuesPage,SelectDimensionsPage,trueLabelDict,DataSelectionHomePage,finalSwitchPage,folderName,secondaryhomepage,processType)
    #master.switch_frame(selectLevelValuesPage,SelectDimensionsPage,trueLabelDict,DataSelectionHomePage,backpage,finalSwitchPage,folderName,secondaryhomepage,processType)
    def __init__(self, master,switchPage,trueLabelDict,backPage,secondaryBackPage,fsp,fName,shp,pt):
        tk.Frame.__init__(self, master)
        
        includeLevelValueList = []
        
        labelWindow = tk.Frame(self)
        labelWindow.pack(side=tk.TOP,padx=10,fill=tk.X,expand=True)
        
        l1 = tk.Label(labelWindow, text='Which specific level values do you want to include in the figure?',pady=10).grid(row=0,column = 0,columnspan=len(trueLabelDict)*6)
        levelValueCheckButtonList = []
        overallCheckButtonVariableList = []
        checkAllButtonList = []
        uncheckAllButtonList = []
        i=0
        maxNumLevelValues = 0
        for levelName in trueLabelDict:
            if len(trueLabelDict[levelName]) > maxNumLevelValues:
                maxNumLevelValues = len(trueLabelDict[levelName])
        """BEGIN TEMP SCROLLBAR CODE"""
        labelWindow1 = tk.Frame(self)
        labelWindow1.pack(side=tk.TOP,padx=10,fill=tk.X,expand=True)
        
        #Make canvas
        w1 = tk.Canvas(labelWindow1, width=1200, height=400,background="white", scrollregion=(0,0,2000,33*maxNumLevelValues))

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
        labelWindow = tk.Frame(w1)
        labelWindow.pack()
        w1.create_window((0,0),window=labelWindow, anchor = tk.NW)
        """END TEMP SCROLLBAR CODE"""
        for levelName in trueLabelDict:
            j=0
            levelCheckButtonList = []
            levelCheckButtonVariableList = []
            levelLabel = tk.Label(labelWindow, text=levelName+':')
            levelLabel.grid(row=1,column = i*6,sticky=tk.N,columnspan=5)
            for levelValue in trueLabelDict[levelName]:
                includeLevelValueBool = tk.BooleanVar()
                cb = tk.Checkbutton(labelWindow, text=levelValue, variable=includeLevelValueBool)
                cb.grid(row=j+4,column=i*6+2,columnspan=2,sticky=tk.W)
                labelWindow.grid_columnconfigure(i*6+3,weight=1)
                cb.select()
                levelCheckButtonList.append(cb)
                levelCheckButtonVariableList.append(includeLevelValueBool)
                j+=1
            
            checkAllButton1 = checkUncheckAllButton(labelWindow,levelCheckButtonList, text='Check All')
            checkAllButton1.configure(command=checkAllButton1.checkAll)
            checkAllButton1.grid(row=2,column=i*6,sticky=tk.N,columnspan=3)
            checkAllButtonList.append(checkAllButton1)
            
            uncheckAllButton1 = checkUncheckAllButton(labelWindow,levelCheckButtonList, text='Uncheck All')
            uncheckAllButton1.configure(command=checkAllButton1.uncheckAll)
            uncheckAllButton1.grid(row=2,column=i*6+3,sticky=tk.N,columnspan=3)
            uncheckAllButtonList.append(checkAllButton1)

            levelValueCheckButtonList.append(levelCheckButtonList)
            overallCheckButtonVariableList.append(levelCheckButtonVariableList)
            i+=1

        def collectInputs():
            for checkButtonVariableList in overallCheckButtonVariableList:
                tempLevelValueList = []
                for checkButtonVariable in checkButtonVariableList:
                    tempLevelValueList.append(checkButtonVariable.get())
                includeLevelValueList.append(tempLevelValueList)
            #master.switch_frame(assignLevelsToParametersPage)
            master.switch_frame(switchPage,includeLevelValueList)
        
        def quitCommand():
            exitBoolean = True
            quit()
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=maxNumLevelValues+4,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backPage,fsp,fName,secondaryBackPage,pt,shp)).grid(row=maxNumLevelValues+4,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quitCommand()).grid(row=maxNumLevelValues+4,column=2)

class assignLevelsToParametersPage(tk.Frame):
    
    def __init__(self, master,temp):
        parameterTypeDict = {
                'categorical':['Color','Order', 'Row', 'Column','None'],
                '1d':['Color','Row','Column','None'],
                '2d':['Marker','Color','Size','Row','Column','X Axis Values','None'],
                '3d':['Row','Column','X Axis Values','Y Axis Values']}
        
        tk.Frame.__init__(self, master)
        global includeLevelValueList
        includeLevelValueList = temp
        global parametersSelected
        parametersSelected = {}
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        l1 = tk.Label(mainWindow, text='Which plotting parameter do you want to assign to each of your figure levels?',pady=10).grid(row=0,column = 0,columnspan = len(figureLevelList))
        rblist = []
        parameterVarList = []
        for figureLevel,figureLevelIndex in zip(figureLevelList,range(len(figureLevelList))):
            v = tk.IntVar()
            temprblist = []
            levelLabel = tk.Label(mainWindow, text=figureLevel+':')
            levelLabel.grid(row=1,column=figureLevelIndex,sticky=tk.NW)
            for plottingParameter,parameterIndex in zip(parameterTypeDict[plotType],range(len(parameterTypeDict[plotType]))):
                rb = tk.Radiobutton(mainWindow, text=plottingParameter,padx = 20, variable=v, value=parameterIndex)
                rb.grid(row=parameterIndex+2,column=figureLevelIndex,sticky=tk.NW)
                temprblist.append(rb)
            rblist.append(temprblist)
            parameterVarList.append(v)
        
        def collectInputs():
            for parameterVar,levelName in zip(parameterVarList,figureLevelList):
                if parameterTypeDict[plotType][parameterVar.get()] not in parametersSelected.keys():
                    parametersSelected[parameterTypeDict[plotType][parameterVar.get()]] = levelName
                else:
                    if not isinstance(parametersSelected[parameterTypeDict[plotType][parameterVar.get()]],list):
                        parametersSelected[parameterTypeDict[plotType][parameterVar.get()]] = [parametersSelected[parameterTypeDict[plotType][parameterVar.get()]]]+[levelName]
                    else:
                        parametersSelected[parameterTypeDict[plotType][parameterVar.get()]].append(levelName)
            master.switch_frame(plotElementsGUIPage)
        
        def quitCommand():
            exitBoolean = True
            quit()
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=len(parameterTypeDict[plotType])+2,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(selectLevelValuesPage,assignLevelsToParametersPage,trueLabelDict,selectLevelsPage,'a','b','c','d','e')).grid(row=len(parameterTypeDict[plotType])+2,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quitCommand()).grid(row=len(parameterTypeDict[plotType])+2,column=2)

def getDefaultAxisTitles():
        xaxistitle = ''
        yaxistitle = ''
        cbartitle = ''
        cytokineDefault = 'Concentration (nM)'
        if plotType == '1d':
            xaxistitle = 'MFI'
            yaxistitle = 'Frequency'
        else:
            if plotType == 'categorical':
                xaxistitle = parametersSelected['Order']
                if dataType == 'cyt':
                    yaxistitle = cytokineDefault
            else:
                xaxistitle = parametersSelected['X Axis Values']
                if dataType == 'cyt':
                    if plotType == '2d':
                        yaxistitle = cytokineDefault
                    else:
                        cbartitle = cytokineDefault
        if xaxistitle == 'Time':
            xaxistitle += ' (hours)'
        return [xaxistitle,yaxistitle,cbartitle]

class plotElementsGUIPage(tk.Frame):
    def __init__(self, master):
        if 'experimentParameters-'+folderName+'-'+expParamDict[dataType]+'.json' in os.listdir('misc'):
            experimentParameters = json.load(open('misc/experimentParameters-'+folderName+'-'+expParamDict[dataType]+'.json','r'))
        else:
            tempDict = {}
            stackedDf = experimentDf.stack()
            for level in stackedDf.index.names:
                levelValues = list(pd.unique(stackedDf.index.get_level_values(level)))
                tempDict[level] = levelValues
            experimentParameters = {}
            experimentParameters['levelLabelDict'] = tempDict
        """ 
        global dataType
        if 'cell' in pickleFileName: 
            dataType = 'cell'
        elif 'cyt' in pickleFileName:
            dataType = 'cyt'
        elif 'prolif' in pickleFileName:
            dataType = 'prolif'
        else:
            dataType = ''
        """
        self.tld = trueLabelDict

        axisDict = {'categorical':['X','Y'],'1d':['Y'],'2d':['X','Y'],'3d':['X','Y','Colorbar']}
        scalingList = ['Linear','Logarithmic','Biexponential']
        axisSharingList = ['col','row','']
        axisTitleDefaults = getDefaultAxisTitles() 
        
        tk.Frame.__init__(self, master)
        
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)

        tk.Label(mainWindow, text='Title: ').grid(row=1,column=0,sticky=tk.W)
        for scaling,scalingIndex in zip(scalingList,range(len(scalingList))):
            tk.Label(mainWindow, text=scaling+' Scaling: ').grid(row=scalingIndex+2,column=0,sticky=tk.W)
        tk.Label(mainWindow, text='Linear Range (Biexponential Scaling): ').grid(row=len(scalingList)+2,column=0,sticky=tk.W)
        tk.Label(mainWindow, text='Convert to numeric: ').grid(row=len(scalingList)+3,column=0,sticky=tk.W)
        tk.Label(mainWindow, text='Share axis across: ').grid(row=len(scalingList)+4,column=0,sticky=tk.W)
        tk.Label(mainWindow, text='Axis limits: ').grid(row=len(scalingList)+5,column=0,sticky=tk.W)

        entryList = []
        scalingVariableList = []
        radioButtonList = []
        checkButtonList = []
        checkButtonVarList = []
        radioButtonList2 = []
        radioButtonVarList2 = []
        linearRangeScalingList = []
        limitEntryList = []
        for axis,axisIndex in zip(axisDict[plotType],range(len(axisDict[plotType]))):
            tk.Label(mainWindow, text=axis+ ' Axis').grid(row=0,column=axisIndex+1)
            
            e1 = tk.Entry(mainWindow)
            e1.grid(row=1,column=axisIndex+1)
            e1.insert(0, axisTitleDefaults[axisIndex])
            entryList.append(e1)
            
            axisRadioButtonList = []
            v = tk.StringVar(value='Linear')
            for scaling,scalingIndex in zip(scalingList,range(len(scalingList))):
                rb = tk.Radiobutton(mainWindow,variable=v,value=scaling)
                rb.grid(row=scalingIndex+2,column=axisIndex+1)
                axisRadioButtonList.append(rb)
            radioButtonList.append(axisRadioButtonList)
            scalingVariableList.append(v)
            
            e2 = tk.Entry(mainWindow)
            e2.grid(row=len(scalingList)+2,column=axisIndex+1)
            linearRangeScalingList.append(e2)
            
            b = tk.BooleanVar(value=False)
            cb = tk.Checkbutton(mainWindow,variable=b)
            cb.grid(row=len(scalingList)+3,column=axisIndex+1)
            checkButtonList.append(cb)
            checkButtonVarList.append(b)
       
            shareWindow = tk.Frame(mainWindow)
            shareWindow.grid(row=len(scalingList)+4,column=axisIndex+1)
            shareString = tk.StringVar(value='None')
            rb2a = tk.Radiobutton(shareWindow,variable=shareString,text='All',value='all')
            rb2b = tk.Radiobutton(shareWindow,variable=shareString,text='Row',value='row')
            rb2c = tk.Radiobutton(shareWindow,variable=shareString,text='Col',value='col')
            rb2d = tk.Radiobutton(shareWindow,variable=shareString,text='None',value='none')
            shareString.set('all')
            rb2a.grid(row=0,column=1)
            rb2b.grid(row=0,column=2)
            rb2c.grid(row=0,column=3)
            rb2d.grid(row=0,column=4)
            radioButtonList2.append([rb2a,rb2b])
            radioButtonVarList2.append(shareString)

            limitWindow = tk.Frame(mainWindow)
            limitWindow.grid(row=len(scalingList)+5,column=axisIndex+1)
            #ll = tk.Label(limitWindow,text='Lower:').grid(row=0,column=0)
            e3 = tk.Entry(limitWindow,width=5)
            e3.grid(row=0,column=1)
            #ul = tk.Label(limitWindow,text='Upper:').grid(row=0,column=2)
            e4 = tk.Entry(limitWindow,width=5)
            e4.grid(row=0,column=3)
            limitEntryList.append([e3,e4])

        def collectInputs(plotAllVar):
            plotOptions = {}
            for axis,axisIndex in zip(axisDict[plotType],range(len(axisDict[plotType]))):
                share = radioButtonVarList2[axisIndex].get()
                if share == 'none':
                    share = False
                plotOptions[axis] = {'axisTitle':entryList[axisIndex].get(),
                        'axisScaling':scalingVariableList[axisIndex].get(),
                        'linThreshold':linearRangeScalingList[axisIndex].get(),
                        'numeric':checkButtonVarList[axisIndex].get(),
                        'share':share,
                        'limit':[limitEntryList[axisIndex][0].get(),limitEntryList[axisIndex][1].get()]}
            
            plotSpecificDict = {}
            if subPlotType == 'kde':
                scaleBool = ipe.getRadiobuttonValues(modeScaleRadiobuttonVarsDict)['scale to mode']
                if scaleBool == 'yes':
                    plotSpecificDict['scaleToMode'] = True
                else:
                    plotSpecificDict['scaleToMode'] = False
                plotSpecificDict['smoothing'] = int(ipe.getSliderValues(smoothingSliderList,['smoothing'])['smoothing'])
            
            useModifiedDf = False
            sName = titleEntry.get()
            subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually = fpl.produceSubsettedDataFrames(experimentDf.stack().to_frame('temp'),fullFigureLevelBooleanList,includeLevelValueList,self.tld)
            fpl.plotFacetedFigures(folderName,plotType,subPlotType,dataType,subsettedDfList,subsettedDfListTitles,figureLevels,levelValuesPlottedIndividually,useModifiedDf,experimentDf,plotOptions,parametersSelected,addDistributionPoints,originalLevelValueOrders=experimentParameters['levelLabelDict'],subfolderName=sName,context=ipe.getRadiobuttonValues(contextRadiobuttonVarsDict)['context'],height=float(heightEntry.get()),aspect=float(widthEntry.get()),titleBool=ipe.getRadiobuttonValues(plotTitleRadiobuttonVarsDict)['plotTitle'],colwrap=int(colWrapEntry.get()),legendBool=ipe.getRadiobuttonValues(legendRadiobuttonVarsDict)['legend'],cmap=cmapEntry.get(),plotAllVar=plotAllVar,titleAdjust=titleAdjustEntry.get(),plotSpecificDict = plotSpecificDict)
            
        titleWindow = tk.Frame(self)
        titleWindow.pack(side=tk.TOP,pady=10)
        tk.Label(titleWindow,text='Enter subfolder for these plots (optional):').grid(row=0,column=0)
        titleEntry = tk.Entry(titleWindow,width=15)
        titleEntry.grid(row=0,column=1)
        
        miscOptionsWindow = tk.Frame(self)
        miscOptionsWindow.pack(side=tk.TOP,pady=10)
        
        contextWindow = tk.Frame(miscOptionsWindow)
        contextWindow.grid(row=0,column=0,sticky=tk.N)
        contextRadiobuttonList,contextRadiobuttonVarsDict = ipe.createParameterSelectionRadiobuttons(contextWindow,['context'],{'context':['notebook','talk','poster']}) 
        
        figureDimensionWindow = tk.Frame(miscOptionsWindow)
        figureDimensionWindow.grid(row=0,column=1,sticky=tk.N)
        tk.Label(figureDimensionWindow,text='figure dimensions').grid(row=0,column=0)
        tk.Label(figureDimensionWindow,text='height:').grid(row=1,column=0)
        tk.Label(figureDimensionWindow,text='width:').grid(row=2,column=0)
        heightEntry = tk.Entry(figureDimensionWindow,width=3)
        if plotType != '1d':
            heightEntry.insert(0, '5')
        else:
            heightEntry.insert(0, '3')
        widthEntry = tk.Entry(figureDimensionWindow,width=3)
        widthEntry.insert(0, '1')
        heightEntry.grid(row=1,column=1)
        widthEntry.grid(row=2,column=1)
        
        plotTitleWindow = tk.Frame(miscOptionsWindow)
        plotTitleWindow.grid(row=0,column=2,sticky=tk.N)
        plotTitleRadiobuttonList,plotTitleRadiobuttonVarsDict = ipe.createParameterSelectionRadiobuttons(plotTitleWindow,['plotTitle'],{'plotTitle':['yes','no']}) 
        
        legendWindow = tk.Frame(miscOptionsWindow)
        legendWindow.grid(row=0,column=3,sticky=tk.N)
        legendRadiobuttonList,legendRadiobuttonVarsDict = ipe.createParameterSelectionRadiobuttons(legendWindow,['legend'],{'legend':['yes','no']}) 
        
        colWrapWindow = tk.Frame(miscOptionsWindow)
        colWrapWindow.grid(row=0,column=4,sticky=tk.N)
        tk.Label(colWrapWindow,text='column wrap:').grid(row=0,column=0)
        colWrapEntry = tk.Entry(colWrapWindow,width=5)
        colWrapEntry.insert(0, '5')
        colWrapEntry.grid(row=1,column=0)
        
        titleAdjustWindow = tk.Frame(miscOptionsWindow)
        titleAdjustWindow.grid(row=0,column=5,sticky=tk.N)
        tk.Label(titleAdjustWindow,text='title location (% of window):').grid(row=0,column=0)
        titleAdjustEntry = tk.Entry(titleAdjustWindow,width=5)
        titleAdjustEntry.insert(0, '')
        titleAdjustEntry.grid(row=1,column=0)
        
        #outlierWindow = tk.Frame(miscOptionsWindow)
        #outlierWindow.grid(row=0,column=6,sticky=tk.N)
        #outlierRadiobuttonList,outlierRadiobuttonVarsDict = ipe.createParameterSelectionRadiobuttons(outlierWindow,['remove outliers'],{'remove outliers':['yes','no']}) 
        #outlierRadiobuttonVarsDict['remove outliers'].set('no')
        
        cmapWindow = tk.Frame(miscOptionsWindow)
        cmapWindow.grid(row=0,column=6,sticky=tk.N)
        tk.Label(cmapWindow,text='Colormap:').grid(row=0,column=0)
        cmapEntry = tk.Entry(cmapWindow,width=10)
        cmapEntry.grid(row=1,column=0)
        
        if subPlotType == 'kde':
            #Scale to mode button
            subPlotSpecificWindow = tk.Frame(self)
            subPlotSpecificWindow.pack(side=tk.TOP,pady=10)
            modeScaleWindow = tk.Frame(subPlotSpecificWindow)
            modeScaleWindow.grid(row=0,column=0,sticky=tk.N)
            modeScaleRadiobuttonList,modeScaleRadiobuttonVarsDict = ipe.createParameterSelectionRadiobuttons(modeScaleWindow,['scale to mode'],{'scale to mode':['yes','no']}) 
            modeScaleRadiobuttonVarsDict['scale to mode'].set('no')
            #Smoothness (or bandwidth) slider
            smoothnessWindow = tk.Frame(subPlotSpecificWindow)
            smoothnessWindow.grid(row=0,column=1,sticky=tk.N)
            smoothingSliderList = ipe.createParameterAdjustmentSliders(smoothnessWindow,['smoothing'],{'smoothing':[1,99,2,27]})
        
        plotButtonWindow = tk.Frame(self)
        plotButtonWindow.pack(side=tk.TOP,pady=10)
        tk.Button(plotButtonWindow, text="Generate First Plot",command=lambda: collectInputs(False)).grid(row=0,column=0)
        tk.Button(plotButtonWindow, text="Generate All Plots",command=lambda: collectInputs(True)).grid(row=0,column=1)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)
        
        def okCommand():
            master.switch_frame(PlotTypePage)
        
        tk.Button(buttonWindow, text="Finish",command=lambda: okCommand()).grid(row=len(scalingList)+4,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(assignLevelsToParametersPage,includeLevelValueList)).grid(row=len(scalingList)+4,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=len(scalingList)+4,column=2)
