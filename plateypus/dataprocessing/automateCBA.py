#! /usr/bin/env python3
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly
import plotly.graph_objs as go
idx = pd.IndexSlice
import os,sys,string
sns.set_context('talk')

#For reading in fcs files in python
import flowio
import flowutils

#For clustering
import hdbscan
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
from sklearn.preprocessing import MinMaxScaler
import json

#For more appropriate umap visualization; might need to install them if using google colab
#import datashader as ds
#import holoviews as hv
#import colorcet as cc
#from datashader import transfer_functions as tf
#import holoviews.operation.datashader as hd

def importFCS(tempChannelDict={},path=''):
    plateNames,plateDfList = [],[]
    subPath = 'inputData/fcsFiles/CBA/'
    for fcsFolder in os.listdir(path+subPath):
        if '.DS' not in fcsFolder and 'Wash' not in fcsFolder:
            #Grab relevant channels
            columnsToKeep,channelsToKeep = [],[]
            for fileName in os.listdir(path+subPath+fcsFolder):
                if '.DS' not in fileName:
                    fd = flowio.FlowData(path+subPath+fcsFolder+'/'+fileName)
                    for channelNum in fd.channels:
                        channelDict = fd.channels[channelNum]
                        if 'PnS' in channelDict.keys():
                            channelName = channelDict['PnS']
                        else:
                            channelName = channelDict['PnN']
                            if channelName in ['FSC-A','SSC-A']:
                                channelName = channelName[:-2]
                            elif channelName in ['PE ( 561 )-A','APC-A','BV421-A','FITC-A']:
                                channelName = tempChannelDict[channelName]
                        if channelName != 'Time':
                            channelsToKeep.append(channelName)
                            columnsToKeep.append(int(channelNum)-1)
                    break
            #Temporary; will start naming things appropriately soon
            if 'Cytokine Identity' in channelsToKeep:
                if 'Multiplex 1' in channelsToKeep:
                    channelsToKeep[channelsToKeep.index('Multiplex 1')] = 'Multiplex 2'
                channelsToKeep[channelsToKeep.index('Cytokine Identity')] = 'Multiplex 1'

            fcsList,fileList,eventList = [],[],[]
            for fileName in os.listdir(path+subPath+fcsFolder):
                if '.DS' not in fileName:
                    fd = flowio.FlowData(path+subPath+fcsFolder+'/'+fileName)

                    events = np.reshape(fd.events, (-1, fd.channel_count))[:,columnsToKeep]
                    fluoro_indices = [x for x in range(events.shape[1]) if x not in [0,1]]
                    xform_events = flowutils.transforms.logicle(events, fluoro_indices,w=1)

                    fcsList.append(xform_events)
                    fileList+=[fileName]*xform_events.shape[0]
                    eventList+=list(range(1,xform_events.shape[0]+1))

            newMI = pd.MultiIndex.from_arrays([fileList,eventList],names=['FileName','Event'])
            rawDf = pd.DataFrame(np.vstack(fcsList),index=newMI,columns=channelsToKeep).query("FSC > 0 and SSC > 0")
            plateDfList.append(rawDf)
            plateNames.append(fcsFolder)

    rawDf = pd.concat(plateDfList,keys=plateNames,names=['Plate'])
    return rawDf

def beadGate(fcsDf,cutoff=0.04,path='',visualize=False,savePlot=False):
    #Move to lognormal space
    fcsDf['FSC'] = np.log10(fcsDf['FSC'])
    fcsDf['SSC'] = np.log10(fcsDf['SSC'])

    #Create min max scaled KDE in each dimension, find value in each dimension where density passes below cutoff, use this to 2D gate the cells
    sizeIntervals = []
    for sizeDimension in ['FSC','SSC']:
        X = fcsDf[sizeDimension]
        space = np.linspace(min(X),max(X),num=1000)
        kde = gaussian_kde(X)
        kdeVals = kde.evaluate(space)
        maxInd = np.argmax(kdeVals)
        minMaxScaledKDE = MinMaxScaler().fit_transform(kdeVals.reshape(-1,1))

        tempDf = pd.DataFrame({sizeDimension:space,'Density':minMaxScaledKDE[:,0]})
        mins = []
        for i in range(maxInd,-1,-1):
            if minMaxScaledKDE[i] < cutoff:
                mins.append(i)
                break
        for i in range(maxInd,len(list(minMaxScaledKDE))):
            if minMaxScaledKDE[i] < cutoff:
                mins.append(i)
                break

        gateInterval = [space[mins[0]],space[mins[1]]]
        sizeIntervals.append(gateInterval)

        if visualize:
            g = sns.relplot(data=tempDf,x=sizeDimension,y='Density',kind='line')
            g.axes.flat[0].axhline(y=cutoff,linestyle=':',color='k')
            g.axes.flat[0].axvline(x=gateInterval[0],linestyle=':',color='k')
            g.axes.flat[0].axvline(x=gateInterval[1],linestyle=':',color='k')
    
    fcsDf['Beads'] = ['Yes' if sizeIntervals[0][0] <= fsc < sizeIntervals[0][1] and sizeIntervals[1][0] <= ssc < sizeIntervals[1][1] else 'No' for fsc,ssc in zip(fcsDf['FSC'],fcsDf['SSC'])]
    """
    if visualize:
        g = facetedSingleCellScatter(data=fcsDf.reset_index(),x='FSC',y='SSC',hue='Beads',hue_order=['Yes','No'],biExpXYScale=False)
        for i in range(2):
            g.axes[0].axvline(x=sizeIntervals[0][i],linestyle=':',color='k')
        for i in range(2):
            g.axes[0].axhline(y=sizeIntervals[1][i],linestyle=':',color='k')
        display(g)
    if savePlot:
        g = facetedSingleCellScatter(data=fcsDf.reset_index(),x='FSC',y='SSC',hue='Beads',hue_order=['Yes','No'],biExpXYScale=False)
        for i in range(2):
            g.axes[0].axvline(x=sizeIntervals[0][i],linestyle=':',color='k')
        for i in range(2):
            g.axes[0].axhline(y=sizeIntervals[1][i],linestyle=':',color='k')
        plt.savefig(path+'plots/automatedCBAProcessingPlots/beadGate.png',bbox_inches='tight')   
    """

    return fcsDf#.query("Beads == 'Yes'").drop('Beads',axis=1)

def cytokineMultiplexGates(fcsDf,multiplexLabelMatrix,multiplexLabels,minClusterSize=10000,minSamples='',path='',visualize=False,savePlot=False):
    
    print(fcsDf.shape)
    if len(multiplexLabels) == 1:
        X = fcsDf.loc[:,multiplexLabels].values.reshape(-1,1)
    else:
        X = fcsDf.loc[:,multiplexLabels]
    cytokineList = [x for x in np.array(multiplexLabelMatrix).flat if x != '']
    numCytokines = len(cytokineList)
    
    if minSamples == '':
        minSamples = int(minClusterSize*0.01)
    
    clusterer = hdbscan.HDBSCAN(min_cluster_size=minClusterSize,min_samples=minSamples)
    clusterer.fit(X.astype(float))
    cluster_labels = list(map(str,clusterer.labels_.tolist()))
    
    fcsDf['Cytokine'] = cluster_labels
    largestClusters = fcsDf['Cytokine'].value_counts().index[:numCytokines].tolist()
    
    centerDf = fcsDf.groupby(['Cytokine'])[multiplexLabels].apply(lambda x: np.mean(x, axis=0)).query("Cytokine == @largestClusters")
    
    if len(multiplexLabels) == 1:
        sortedDf = centerDf.sort_values(by=multiplexLabels[0])
        splitCytokines = sortedDf.index.tolist()
    elif len(multiplexLabels) == 2:
        sortedDf = centerDf.sort_values(by=multiplexLabels[1],ascending=False)
        splitMatrices2,splitCytokines = [],[]
        splitIndexStart = 0
        for splitMatrix in np.split(multiplexLabelMatrix,multiplexLabelMatrix.shape[0],axis=0):
            numSplitCytokines = splitMatrix[splitMatrix != ''].shape[1]
            splitDf2 = sortedDf.iloc[splitIndexStart:splitIndexStart+numSplitCytokines,:]
            splitCytokines+=splitDf2.sort_values(by=multiplexLabels[0]).index.tolist()
            splitIndexStart+=numSplitCytokines
    
    cytokineDict = {**{key:value for (key,value) in zip(splitCytokines,cytokineList)},**{'-1':'Noise'}}
    
    cytokineList2 = [str(cytokineDict[x]) if x in list(cytokineDict.keys()) else x for x in fcsDf['Cytokine']]
    fcsDf = fcsDf.drop('Cytokine',axis=1).assign(Cytokine=cytokineList2).set_index('Cytokine', append=True).swaplevel(-1,-2)
    """
    if visualize:
        if len(multiplexLabels) == 1:
            g = sns.displot(data=fcsDf.reset_index(),x=multiplexLabels[0],hue='Cytokine',kind='kde',palette=cc.glasbey[:numCytokines+1],hue_order=['Noise']+cytokineList)
        if len(multiplexLabels) == 2:
            #g = sns.relplot(data=fcsDf,x=multiplexLabels[0],y=multiplexLabels[1],hue='Cytokine',s=1,alpha=0.7,palette=cc.glasbey[:numCytokines+1],hue_order=['Noise']+cytokineList)   
            g = facetedSingleCellScatter(data=fcsDf.fillna(value=0).astype(float).reset_index(),x=multiplexLabels[0],y=multiplexLabels[1],hue='Cytokine',biExpXYScale=False,palette=cc.glasbey[:numCytokines+1],hue_order=['Noise']+cytokineList)
            display(g)
    if savePlot:
        if len(multiplexLabels) == 1:
            g = sns.displot(data=fcsDf.reset_index(),x=multiplexLabels[0],hue='Cytokine',kind='kde',palette=cc.glasbey[:numCytokines+1],hue_order=['Noise']+cytokineList)
        if len(multiplexLabels) == 2:
            g = facetedSingleCellScatter(data=fcsDf.fillna(value=0).astype(float).reset_index(),x=multiplexLabels[0],y=multiplexLabels[1],hue='Cytokine',biExpXYScale=False,palette=cc.glasbey[:numCytokines+1],hue_order=['Noise']+cytokineList)
        plt.savefig(path+'plots/automatedCBAProcessingPlots/multiplexGate.png',bbox_inches='tight')   
    """

    return fcsDf.query("Cytokine != 'Noise'")

def plateBarcodeGate(fcsDf,barcodes,barcodingDict='',path='',visualize=False,savePlot=False):
    if len(barcodes) > 0:
        #Beads/IFNg/Q1: Barcode 1- , Barcode 2+ | Geometric Mean (Cytokine Level)
        barcodeClusterDfList = []
        barcodeList = []
        barcodeName = ','.join(barcodes)
        for cluster in fcsDf.index.unique('Cytokine'):
            tempDf = fcsDf.query("Cytokine == @cluster")
            clusteringMatrix = tempDf.loc[:,barcodes].fillna(value=0).values    
            trueMins = []
            for col in range(clusteringMatrix.shape[1]):
                X = clusteringMatrix[:,col]
                space = np.linspace(min(X),max(X),num=1000)
                kde = gaussian_kde(X)
                kdeVals = kde.evaluate(space)
                maxInd = argrelextrema(np.nan_to_num(kdeVals,nan=0), np.greater)
                minInd = argrelextrema(np.nan_to_num(kdeVals,nan=0), np.less)

                # get the actual values using these indices
                relMaxes = space[maxInd]  # array([5, 3, 6])
                maxVals = kde.evaluate(relMaxes).tolist()
                trueMaxes = [relMaxes[maxVals.index(x)] for x in sorted(maxVals,reverse=True)]
                trueMaxVals = [maxVals[maxVals.index(x)] for x in sorted(maxVals,reverse=True)]
                peakRatio = trueMaxVals[0]/trueMaxVals[1]
                peakRatioCutoff = 5
                if len(trueMaxes) >= 2 and peakRatio < peakRatioCutoff:
                    trueMaxes = sorted(trueMaxes[:2])
                    relMins = space[minInd]  # array([5, 3, 6])
                    trueMin = 0
                    for relMin in relMins:
                        if trueMaxes[0] < relMin < trueMaxes[1]: 
                            trueMin = relMin
                            break
                    trueMins.append(trueMin)
                else:
                    trueMins.append(space[500])

            rawBarcodingMatrix = np.empty([tempDf.shape[0],len(barcodes)],dtype=str)
            for i,barcode in enumerate(barcodes):
                rawBarcodingMatrix[:,i] = ['+' if x == 1 else '-' for x in (tempDf.loc[:,barcode].values > trueMins[i]).astype(int)]

            barcodeList+=[','.join(row) for row in rawBarcodingMatrix]
            barcodeClusterDfList.append(tempDf)

        fcsDf = pd.concat(barcodeClusterDfList)
        fcsDf = fcsDf.assign(Barcodes=barcodeList).set_index('Barcodes', append=True)
        
        """
        if visualize:
            if len(barcodes) == 1:
                g = sns.displot(data=fcsDf.astype(float).reset_index(),x=barcodes[0],hue='Barcodes',col='Cytokine',kind='kde',col_wrap=4)    
            elif len(barcodes) == 2:
                g = facetedSingleCellScatter(data=fcsDf.astype(float).reset_index(),x=barcodes[0],y=barcodes[1],biExpXYScale=False,hue='Barcodes',col='Cytokine',col_wrap=4)
                display(g)
        if savePlot:
            if len(barcodes) == 1:
                g = sns.displot(data=fcsDf.astype(float).reset_index(),x=barcodes[0],hue='Barcodes',col='Cytokine',kind='kde',col_wrap=4)    
            elif len(barcodes) == 2:
                g = facetedSingleCellScatter(data=fcsDf.astype(float).reset_index(),x=barcodes[0],y=barcodes[1],biExpXYScale=False,hue='Barcodes',col='Cytokine',col_wrap=4)
            plt.savefig(path+'plots/automatedCBAProcessingPlots/barcodeGate.png',bbox_inches='tight')   
        """
        #Only select used barcode combinations
        if barcodingDict != '':
            newBarcodesToKeep = []
            for key in barcodingDict:
                val = ','.join(barcodingDict[key])
                if val not in newBarcodesToKeep:
                    newBarcodesToKeep.append(val)
            finalBarcodesToKeep = []
            for barcodeToKeep in newBarcodesToKeep:
                b = barcodeToKeep.split(',')
                finalBarcodesToKeep.append(','.join([x[-1] for x in b]))
            fcsDf = fcsDf.query("Barcodes == @finalBarcodesToKeep")

        fcsDf.index.names = [barcodeName if x == 'Barcodes' else x for x in fcsDf.index.names]
        fcsDf = fcsDf.swaplevel(-1,-2)
    
    return fcsDf

#Beads/IFNg/Q1: Barcode 1- , Barcode 2+ | Geometric Mean (Cytokine Level)
#Beads/IFNg/Barcode 1- | Geometric Mean (Cytokine Level)
#Beads/IFNg | Geometric Mean (Cytokine Level)
def produceBulkCSVFiles(completeFCSDf,path=''):
    for plateName in completeFCSDf.index.unique('Plate'):
        if 'Calibration' in plateName and '-' in plateName:
            beadsName = 'Beads/'+plateName.split('-')[1]
            #beadsName = 'Beads'
        else:
            beadsName = 'Beads'
        fcsDf = completeFCSDf.xs([plateName],level=['Plate'])
        if len(fcsDf.index.names) == 4:
            unstackingLevels = list(fcsDf.index.names)[1:]
        else:
            unstackingLevels = [list(fcsDf.index.names)[1]]       
        fcsDf = fcsDf.groupby(['FileName']+unstackingLevels).mean().loc[:,'Cytokine Level'].unstack(unstackingLevels)
        fcsDf.iloc[:,:] = flowutils.transforms.logicle_inverse(fcsDf.values, list(range(fcsDf.shape[1])),w=1)
        suffix = ' | Geometric Mean (Cytokine Level)'
        #Barcoding
        if len(fcsDf.index.names) == 4:
            cytokineIndex = [x for x in range(fcsDf.columns.names) if x == 'Cytokine'][0]
            barcodeIndex = [x for x in range(fcsDf.columns.names) if x != 'Cytokine'][0]
            newColumnList = []
            for col in range(len(fcsDf.columns)):
                name = list(fcsDf.iloc[:,col].name)
                cytokine = name[cytokineIndex]
                barcode = name[barcodeIndex]
                barcodeList = []
                for i,b in enumerate(barcode.split(',')):
                    barcodeList.append('Barcode '+str(i+1)+b)
                population = '/'.join([beadsName,cytokine]+barcodeList)
                newColumnList.append(population+suffix)
            fcsDf.columns = newColumnList
        #No barcoding
        else:
            prefix = beadsName+'/'
            fcsDf.columns = [prefix+x+suffix for x in fcsDf.columns]
        
        specimenOrderingDict = {}
        for i,letter in enumerate(string.ascii_uppercase[:16]):
            for number in range(24):
                trueNumber = i*24 + number
                specimenOrderingDict[letter+str(number+1).zfill(2)] = trueNumber
        csvOrder = [specimenOrderingDict[x.split('.')[0].split('_')[-2]] for x in fcsDf.index.get_level_values('FileName').tolist()]
        finalUnstackedCSV = fcsDf.assign(Order=csvOrder).set_index('Order', append=True).sort_index(level='Order').droplevel('Order')
        finalUnstackedCSV.columns.name = ''
        finalUnstackedCSV.index.name = ''
        
        suffixDf = pd.DataFrame(np.ones([2,finalUnstackedCSV.shape[1]]),index=['Mean','SD'],columns=finalUnstackedCSV.columns)
        finalUnstackedCSV = pd.concat([finalUnstackedCSV,suffixDf]).round().fillna(value=1).astype(int)      
        finalUnstackedCSV[''] = np.nan

        name = plateName+'_cyt.csv'
        finalUnstackedCSV.to_csv(path+'inputData/bulkCSVFiles/'+name)

def automatedCBAProcessing(multiplexes,barcodes,cytokineList,tempChannelDict='',barcodingDict='',cutoff=0.04,minClusterSize = '',path='',name='',visualize=False,savePlot=False,saveDf=False):
    rawDf = importFCS(tempChannelDict=tempChannelDict,path=path)
    print('FCS Files Imported')
    beadGatedDf = beadGate(rawDf,cutoff=cutoff,path=path,visualize=visualize,savePlot=savePlot)
    print('Bead Gate Completed')
    #Automatic HDBScan cluster size collection
    if minClusterSize == '':
        minClusterSize = int(beadGatedDf.shape[0]/len(cytokineList)/1.5)
    cytokineGatedDf = cytokineMultiplexGates(beadGatedDf,np.matrix(cytokineList),multiplexes,minClusterSize=minClusterSize,path=path,visualize=visualize,savePlot=savePlot)
    print('Cytokine Gate Completed')
    barcodeGatedDf = plateBarcodeGate(cytokineGatedDf,barcodes,barcodingDict=barcodingDict,path=path,visualize=visualize,savePlot=savePlot)
    print('Barcoding Gate Completed')
#    sys.exit(0)
    produceBulkCSVFiles(barcodeGatedDf,path=path)
    if saveDf:
        barcodeGatedDf.to_pickle(path+'outputData/pickleFiles/singleBeadDf-'+name+'.pkl')
