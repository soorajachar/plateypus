#!/usr/bin/env python3 
import pandas as pd
from collections import Counter
import string
import sys

def parseCellCSVHeaders(columns,panelData=[]):    
    #Determine marker and statistic first; simply add the population name as is for now
    newMultiIndexList = []
    for column in columns:
        if 'Unnamed' not in column:
            populationNameVsStatisticSplit = column.split(' | ')
            fullPopulationName = populationNameVsStatisticSplit[0]
            #Statistics can be performed on the whole cell population, in which case the cellType is allEvents
            #MFI and CV need to be specified in terms of a laser channel; count and percent positive do not
            if '/' in fullPopulationName:
                populationDivisionIndices = [i for i,c in enumerate(fullPopulationName) if c=='/']
                cellType = fullPopulationName
                #MFI or CV
                #print(populationNameVsStatisticSplit[1])
                if 'Mean' in populationNameVsStatisticSplit[1] or 'Median' in populationNameVsStatisticSplit[1] or 'CV' in populationNameVsStatisticSplit[1]:
                    cellType = fullPopulationName
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
                #% of parent and count
                else:
                    marker = 'NotApplicable'
                    if 'Freq' in populationNameVsStatisticSplit[1]:
                        statistic = '%'
                    else:
                        statistic = 'Count'
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
    
    #Now determine population name
    commonBranches = []
    allFullPopNames = []
    for multiIndexList in newMultiIndexList:
        fullPopulationName = multiIndexList[0]
        commonBranches.append(fullPopulationName.split('/'))
        allFullPopNames.append(fullPopulationName)
    
    #Only one unique sequence of gates; simply take last gate name
    if len(set(allFullPopNames)) == 1:
        populationNames = [allFullPopNames[0].split('/')[-1]]*len(newMultiIndexList)
    #Otherwise use last common gate name as root of all other gating trees
    else:
        minGatingTreeLen = len(min(commonBranches, key=len))
        commonGates,commonGateIndices = [],[]
        for gateIndex in range(minGatingTreeLen):
            currentLevel = []
            for cb in commonBranches:
                currentLevel.append(cb[gateIndex])
            if len(set(currentLevel)) == 1:
                commonGates.append(cb[gateIndex])
                commonGateIndices.append(gateIndex)
            else:
                break
        #If there is at least one common population name
        if len(commonGateIndices) >= 1:
            lastCommonGateIndex = commonGateIndices[-1]
            populationNames = ['/'.join(x.split('/')[lastCommonGateIndex:]) for x in allFullPopNames]
        else:
            populationNames = allFullPopNames.copy()

    #Add cropped population names to multiIndex
    for i,populationName in enumerate(populationNames):
        newMultiIndexList[i] = [populationName]+newMultiIndexList[i][1:]

    trueMIList = [';'.join(x) for x in newMultiIndexList]
    trueUniqueMIList = list(set(trueMIList))
    cnt = Counter(trueMIList)
    nonUniqueItems = [k for k, v in cnt.items() if v > 1]
    
    if len(nonUniqueItems) > 0:
        uppercase = string.ascii_uppercase
        printingNonUnique = [x.split(';') for x in nonUniqueItems]
        print('These column headers are duplicated: please remove them from the table editor, re-export, and try reprocessing:')
        for p in printingNonUnique:
            indexVal = trueMIList.index(';'.join(p))+2
            n = indexVal
            result = ''
            while n > 0:
                index = (n-1)%26
                result = uppercase[int(index)]+result
                n = int((n-1)/26)
            print('Column '+str(indexVal)+' ('+result+'): '+str(p))

    return newMultiIndexList

def parseCellCSVHeadersOld(columns,panelData=[]):

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
    print(newMultiIndexList)
    return newMultiIndexList

def createCellDataFrame(folderName,finalDataFrame,concUnitPrefix):
    return finalData
