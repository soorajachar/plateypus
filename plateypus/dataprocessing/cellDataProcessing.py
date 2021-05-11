#!/usr/bin/env python3 
import pickle
import pandas as pd
import tkinter as tk
from plateypus.dataprocessing.miscFunctions import setMaxWidth

FACS_Detector_Names = ['BV421-A','BV510-A','BV605-A','BV605-A','BV650-A','BV711-A','BV786-A','BUV396-A','DAPI-A','BUV737-A','APC-A','Alexa Fluor 700-A','APC-Cy7-A','FITC-A','PerCP-Cy5-5-A','PE ( 561 )-A','PE-CF594-A','PE-Cy5-A','PE-Cy7-A']
CyTOF_Detector_Names = []
Full_Detector_Dict = {'FACS':FACS_Detector_Names,'CyTOF':CyTOF_Detector_Names}

class MarkerNumberPage(tk.Frame):
    def __init__(self, master,folderName,expNum,ex_data,shp,bPage):
        tk.Frame.__init__(self, master)
        
        global secondaryhomepage
        secondaryhomepage = shp

        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10) 

        l1 = tk.Label(mainWindow,text='Cytometer Type: ')
        l1.grid(row=0,column=0,sticky=tk.W)
        v = tk.StringVar(value='FACS')
        rb1a = tk.Radiobutton(mainWindow,text='FACS',variable=v,value='FACS')
        rb1b = tk.Radiobutton(mainWindow,text='CyTOF',variable=v,value='CyTOF')
        rb1a.grid(row=1,column=0,sticky=tk.W)
        rb1b.grid(row=2,column=0,sticky=tk.W)

        l2 = tk.Label(mainWindow,text='Number of Markers: ')
        l2.grid(row=0,column=1,sticky=tk.W)
        t = tk.Entry(mainWindow)
        t.grid(row=1,column=1,sticky=tk.W)

        def collectInputs():
            numMarkers = int(t.get())
            cytometerType = v.get()
            master.switch_frame(AntibodyPanelPage,numMarkers,cytometerType,folderName,expNum,ex_data,bPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

class AntibodyPanelPage(tk.Frame):
    def __init__(self, master,numMarkers,cytometerType,folderName,expNum,ex_data):
        tk.Frame.__init__(self, master)
        
        experimentNameWindow = tk.Frame(self)
        experimentNameWindow.pack(side=tk.TOP,padx=10,pady=10)
        experimentNameLabel = tk.Label(experimentNameWindow,text=folderName+':').pack()

        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10) 
        
        l1 = tk.Label(mainWindow,text='Marker: ').grid(row=0,column=1,sticky=tk.W)
        l2 = tk.Label(mainWindow,text='Channel: ').grid(row=0,column=0,sticky=tk.W)
        l3 = tk.Label(mainWindow,text='T Cell Gate: ').grid(row=0,column=2,sticky=tk.W)
        l4 = tk.Label(mainWindow,text='APC Gate: ').grid(row=0,column=3,sticky=tk.W)

        markerEntryList = []
        dropMenuList = []
        dropVarList = []
        cblist = []
        cbvarlist = []

        for marker in range(numMarkers):
            v = tk.StringVar()
            t = tk.Entry(mainWindow)
            t.grid(row=marker+1,column=1)
            markerEntryList.append(t)
            dropVar = tk.StringVar()
            dropMenu = tk.OptionMenu(mainWindow,dropVar,*Full_Detector_Dict[cytometerType])
            dropMenu.grid(row=marker+1,column=0,sticky=tk.EW)
            setMaxWidth(Full_Detector_Dict[cytometerType],dropMenu)
            dropVarList.append(dropVar)
            dropMenuList.append(dropMenu)
            cbvar1 = tk.BooleanVar()
            cb1 = tk.Checkbutton(mainWindow,variable=cbvar1)
            cbvar2 = tk.BooleanVar()
            cb2 = tk.Checkbutton(mainWindow,variable=cbvar2)
            cb1.grid(row=marker+1,column=2,sticky=tk.N)
            cb2.grid(row=marker+1,column=3,sticky=tk.N)
            cblist.append([cb1,cb2])
            cbvarlist.append([cbvar1,cbvar2])
        
        def collectInputs():
            finalMarkerList = []
            gatingDict = {}
            finalChannelList = []
            for marker in range(numMarkers):
                if cbvarlist[marker][0].get():
                    gatingDict['TCell_Gate'] = markerEntryList[marker].get()
                elif cbvarlist[marker][1].get():
                    gatingDict['APC_Gate'] = markerEntryList[marker].get()
                finalMarkerList.append(markerEntryList[marker].get())
                finalChannelList.append(dropVarList[marker].get())
            dfDict = {'Marker':finalMarkerList,'FCSDetectorName':finalChannelList}
            antibodyPanelDataFrame = pd.DataFrame(dfDict)
            antibodyPanelDataFrame.to_csv('inputFiles/antibodyPanel-'+folderName+'.csv', encoding='utf-8', index=False)
            with open('misc/gateVals.pkl','wb') as f:
                pickle.dump(gatingDict,f)
            master.switch_frame(secondaryhomepage,folderName,expNum,ex_data)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,pady=10)

        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs()).grid(row=5,column=0)
        tk.Button(buttonWindow, text="Quit",command=quit).grid(row=5,column=2)

def parseCellCSVHeaders(columns,panelData=[]):

#,SingleCells/CD45Neg/TumorCells | Count,SingleCells/CD45Neg/TumorCells | Geometric Mean (Comp-APC-A),SingleCells/CD45Neg/TumorCells | Geometric Mean (Comp-BV650-A),SingleCells/CD45Neg/TumorCells | Geometric Mean (Comp-PE ( 561 )-A),SingleCells/CD45Pos/TCells | Count,SingleCells/CD45Pos/TCells | Geometric Mean (Comp-APC-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-APC-Cy7-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-BUV737-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-BV421-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-BV510-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-BV650-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-BV711-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-BV786-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-FITC-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-PE-CF594-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-PE-Cy7-A),SingleCells/CD45Pos/TCells | Geometric Mean (Comp-PerCP-Cy5-5-A),SingleCells/CD45Pos/TCells/CD25+ | Freq. of Parent (%),SingleCells/CD45Pos/TCells/CD25+ | Count,SingleCells/CD45Pos/TCells/CD27+ | Freq. of Parent (%),SingleCells/CD45Pos/TCells/CD27+ | Count,SingleCells/CD45Pos/TCells/CD39+ | Freq. of Parent (%),SingleCells/CD45Pos/TCells/CD39+ | Count,SingleCells/CD45Pos/TCells/CD44+ | Freq. of Parent (%),SingleCells/CD45Pos/TCells/CD44+ | Count,SingleCells/CD45Pos/TCells/CD54+ | Freq. of Parent (%),SingleCells/CD45Pos/TCells/CD54+ | Count,SingleCells/CD45Pos/TCells/CD69+ | Freq. of Parent (%),SingleCells/CD45Pos/TCells/CD69+ | Count,SingleCells/CD45Pos/TCells/PD1+ | Freq. of Parent (%),SingleCells/CD45Pos/TCells/PD1+ | Count,
    #,Cells/Single Cells/APCs | Geometric Mean (Comp-BV605-A),Cells/Single Cells/APCs | Geometric Mean (Comp-FITC-A),Cells/Single Cells/APCs | Geometric Mean (Comp-PE ( 561 )-A),Cells/Single Cells/APCs | Count,Cells/Single Cells/APCs/CD86+ | Freq. of Parent (%),Cells/Single Cells/APCs/H2Kb+ | Freq. of Parent (%),Cells/Single Cells/APCs/PDL1+ | Freq. of Parent (%),Cells/Single Cells/TCells | Count,Cells/Single Cells/TCells | Geometric Mean (Comp-BUV737-A),Cells/Single Cells/TCells | Geometric Mean (Comp-BV605-A),Cells/Single Cells/TCells | Geometric Mean (Comp-PE-CF594-A),Cells/Single Cells/TCells | Geometric Mean (Comp-PE-Cy7-A),Cells/Single Cells/TCells | Geometric Mean (Comp-PerCP-Cy5-5-A),Cells/Single Cells/TCells/CD27+ | Freq. of Parent (%),Cells/Single Cells/TCells/CD54+ | Freq. of Parent (%),Cells/Single Cells/TCells/CD69+ | Freq. of Parent (%),Cells/Single Cells/TCells/PDL1+ | Freq. of Parent (%),
    newMultiIndexList = []
    for column in columns:
        if 'Unnamed' not in column:
            populationNameVsStatisticSplit = column.split(' | ')
            fullPopulationName = populationNameVsStatisticSplit[0]
            #Statistics can be performed on the whole cell population, in which case the cellType is allEvents
            #MFI and CV need to be specified in terms of a laser channel; count and percent positive do not
            if '/' in fullPopulationName:
                populationDivisionIndices = [i for i,c in enumerate(fullPopulationName) if c=='/']
                #MFI or CV
                #print(populationNameVsStatisticSplit[1])
                if 'Mean' in populationNameVsStatisticSplit[1] or 'Median' in populationNameVsStatisticSplit[1] or 'CV' in populationNameVsStatisticSplit[1]:
                    cellType = fullPopulationName
                    #cellType = fullPopulationName[populationDivisionIndices[-1]+1:]
                    if 'Comp-' in populationNameVsStatisticSplit[1]:
                        statisticVsChannelSplit = populationNameVsStatisticSplit[1].split(' (Comp-')
                    else:
                        statisticVsChannelSplit = populationNameVsStatisticSplit[1].split(' (')
                    statistic = statisticVsChannelSplit[0]
                    #print(statistic)
                    if 'Mean' in statistic:
                        statisticName = 'MFI'
                    elif 'Median' in statistic:
                        statisticName = 'MedianFI'
                    else:
                        statisticName = 'CV'
                    
                    #Can have IRF4+,TCells,CD45RB+CD25-; first case: cellType = first case:cellType=previousPop,marker=currentPop[:-1],stat=Positive/Negative MFI; 
                    #second case: cellType = currentPop,marker=NotApplicable,statistic=MFI; third case: cellType = currentPop,marker=NotApplicable,statistic=MFI
                    numPositive = fullPopulationName[populationDivisionIndices[-1]+1:].count('+')
                    numNegative = fullPopulationName[populationDivisionIndices[-1]+1:].count('-')
                    channel = statisticVsChannelSplit[1][:-1]
                    if '-A' in channel:
                        if '::' not in statisticVsChannelSplit[1]:
                            panelIndex = list(panelData['FCSDetectorName']).index(channel)
                            marker = panelData['Marker'][panelIndex]
                        else:
                            marker = statisticVsChannelSplit[1].split(' :: ')[-1][:-1] 
                    else:
                        if '::' not in statisticVsChannelSplit[1]:
                            marker = channel 
                        else:
                            marker = statisticVsChannelSplit[1].split(' :: ')[-1][:-1] 
                    #if numPositive+numNegative == 0, just use the last gate as the population name:
                    if numPositive+numNegative == 0:
                        cellType = fullPopulationName[populationDivisionIndices[-1]+1:]
                        statistic = statisticName
                    #If more than 1 +/- signs, include full population gating hiearchy as cell population name
                    elif numPositive+numNegative > 1:
                        cellType = fullPopulationName
                        statistic = statisticName
                    #If just a single + or - sign, extract marker name out of last gate
                    else:
                        positivePop = fullPopulationName[populationDivisionIndices[-1]+1:-1]
                        if positivePop == marker:
                            if len(populationDivisionIndices) > 1:
                                cellType = fullPopulationName[populationDivisionIndices[-2]+1:populationDivisionIndices[-1]]
                            else:
                                cellType = fullPopulationName[:populationDivisionIndices[-1]]
                        else:
                            cellType = fullPopulationName
                        if numNegative == 0: 
                            if cellType.split('/')[-1][:-1] == marker:
                                statistic = 'Positive '+statisticName
                            else:
                                statistic = statisticName 
                        else:
                            if cellType.split('/')[-1][:-1] == marker:
                                statistic = 'Negative '+statisticName
                            else:
                                statistic = statisticName
                #% of parent and count
                else:
                    numPositive = fullPopulationName[populationDivisionIndices[-1]+1:].count('+')
                    numNegative = fullPopulationName[populationDivisionIndices[-1]+1:].count('-')
                    if 'Freq' in populationNameVsStatisticSplit[1]:
                        if numPositive+numNegative == 0:
                            marker = 'NotApplicable'
                            statistic ='%'
                            #cellType = fullPopulationName[populationDivisionIndices[-1]+1:]
                            cellType = fullPopulationName
                        elif numPositive+numNegative > 1:
                            marker = 'NotApplicable'
                            statistic = '%'
                            cellType = fullPopulationName
                        else:
                            if numNegative == 0:
                                statistic = '% Positive'
                            else:
                                statistic = '% Negative'
                            marker = fullPopulationName[populationDivisionIndices[-1]+1:len(fullPopulationName)-1]
                            positivePop = fullPopulationName[populationDivisionIndices[-1]+1:-1]
                            if positivePop == marker:
                                if len(populationDivisionIndices) > 1:
                                    cellType = fullPopulationName[populationDivisionIndices[-2]+1:populationDivisionIndices[-1]]
                                else:
                                    cellType = fullPopulationName[:populationDivisionIndices[-1]]
                            else:
                                cellType = fullPopulationName
                    else:
                        marker = 'NotApplicable'
                        statistic = 'Count'
                        if numPositive+numNegative == 0:
                            cellType = fullPopulationName[populationDivisionIndices[-1]+1:] 
                        else:
                            cellType = fullPopulationName
            else:
                #Statistics can be performed on the whole cell population, in which case the cellType is allEvents
                #DAPI+ | Freq. of Parent (%)
                #Positive cell percentage statistics do not have channel names, so treat differently
                #MFI or CV
                if len(populationNameVsStatisticSplit) == 1:
                    cellType = 'allEvents'
                    populationNameVsStatisticSplit = [' ',populationNameVsStatisticSplit[0]]
                else:
                    cellType = populationNameVsStatisticSplit[0]
                if 'Mean' in populationNameVsStatisticSplit[1] or 'Median' in populationNameVsStatisticSplit[1] or 'CV' in populationNameVsStatisticSplit[1]:
                    if 'Comp-' in populationNameVsStatisticSplit[1]:
                        statisticVsChannelSplit = populationNameVsStatisticSplit[1].split(' (Comp-')
                    else:
                        statisticVsChannelSplit = populationNameVsStatisticSplit[1].split(' (')
                    statistic = statisticVsChannelSplit[0]
                    if 'Mean' in statistic:
                        statistic = 'MFI'
                    elif 'Median' in statistic:
                        statistic = 'MedianFI'
                    else:
                        statistic = 'CV'
                    channel = statisticVsChannelSplit[1][:-1]
                    marker = channel
                    #panelIndex = list(panelData['FCSDetectorName']).index(channel)
                    #marker = panelData['Marker'][panelIndex]
                #% of parent and count
                else:
                    marker = 'NotApplicable'
                    if 'Freq' in populationNameVsStatisticSplit[1]:
                        statistic = '% Positive'
                    else:
                        statistic = 'Count'
            newMultiIndexList.append([cellType,marker,statistic])
    
    commonBranchesIndices = []
    commonBranches = []
    for multiIndexList in newMultiIndexList:
        fullPopulationName = multiIndexList[0]
        populationDivisionIndices = [i for i,c in enumerate(fullPopulationName) if c=='/']
        commonBranches.append(fullPopulationName.split('/'))
        commonBranchesIndices.append(populationDivisionIndices)
    commonIndex = 0
    masterBranchList = []
    for branchlist in commonBranches:
        for branch in branchlist:
            masterBranchList.append(branch)
    uniqueBranches = list(pd.unique(masterBranchList))
    commonToAllStatistics = []
    for uniqueBranch in uniqueBranches:
        isCommonToAllStatistics = True
        for branchlist in commonBranches:
            if uniqueBranch not in branchlist:
                isCommonToAllStatistics = False
        if isCommonToAllStatistics:
            commonToAllStatistics.append(uniqueBranch)
    uncommonStatisticsList = []
    for commonBranch in commonBranches:
        uncommonStatistics = []
        for branch in commonBranch:
            if branch not in commonToAllStatistics:
                uncommonStatistics.append(branch)
        uncommonStatisticsList.append('/'.join(uncommonStatistics))
    for uncommonStatistic,i in zip(uncommonStatisticsList,range(len(uncommonStatisticsList))):
        newMultiIndexList[i] = [uncommonStatistic]+newMultiIndexList[i][1:]
    return newMultiIndexList

def createCellDataFrame(folderName,finalDataFrame,concUnitPrefix):
    return finalData
