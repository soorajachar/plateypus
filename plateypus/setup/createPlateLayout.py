#!/usr/bin/env python3 
import pickle,math,matplotlib,sys,os,string,subprocess
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import seaborn as sns
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RectangleSelector
from plateypus.dataprocessing.miscFunctions import rainbow_text

idx = pd.IndexSlice
dataTypeObservableRangeDict = {'cyt':1,'cell':3,'prolif':1}
realDataTypeNameDict = {'cyt':'Supernatant','cell':'Surface/Intracellular Marker','prolif':'Proliferation'}
plateRowLetters = string.ascii_uppercase[:16]
plateColumnNumbers = list(range(1,25))

colwrap = 6 
#For macbook by itself
#figLengthScaling = 0.5
#For work monitor
figLengthScaling = 1
#For home monitor
#figLengthScaling = 0.75

#https://stackoverflow.com/questions/16856788/slice-2d-array-into-smaller-2d-arrays
def blockshaped(arr, nrows, ncols):
    """
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size

    If arr is a 2D array, the returned array should look like n subblocks with
    each subblock preserving the "physical" layout of arr.
    """
    h, w = arr.shape
    assert h % nrows == 0, "{} rows is not evenly divisble by {}".format(h, nrows)
    assert w % ncols == 0, "{} cols is not evenly divisble by {}".format(w, ncols)
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def returnBaseLayout(plateDimensions,conditionPlateRows,timepointPlateColumns):
    if conditionPlateRows == 1 and timepointPlateColumns > colwrap:
        colwrapBool = True
    else:
        colwrapBool = False

    if colwrapBool:
        visualNumRowPlates = math.ceil(timepointPlateColumns/colwrap)
        visualNumColumnPlates = colwrap
    else:
        visualNumRowPlates = conditionPlateRows
        visualNumColumnPlates = timepointPlateColumns

    baseLayoutX = list(range(1,plateDimensions[1]*visualNumColumnPlates+1))
    baseLayoutY = list(range(plateDimensions[0]*visualNumRowPlates,0,-1))
    
    if colwrapBool:
        rowLimits = list(range((visualNumRowPlates-1)*plateDimensions[0],len(baseLayoutY)))
        plateOverhang = timepointPlateColumns - ((visualNumRowPlates-1)*colwrap)
        colLimits = list(range((plateOverhang*plateDimensions[1]),(visualNumColumnPlates*plateDimensions[1])))
    else:
        rowLimits = []
        colLimits = []

    #Slightly shift points over every time plate changes
    hshift = np.repeat(list(range(visualNumColumnPlates)),plateDimensions[1])
    vshift = np.repeat(list(range(visualNumRowPlates)),plateDimensions[0])

    baseLayoutX = np.add(np.array(baseLayoutX),np.array(hshift))
    baseLayoutY = np.add(np.array(baseLayoutY),np.array(vshift[::-1]))

    pointList = []
    infoList = []
    trueCol = -1

    for row,y in enumerate(baseLayoutY):
        if row % plateDimensions[0] == 0:
            trueCol+=1
        for col,x in enumerate(baseLayoutX):
            rowLetter = plateRowLetters[row%plateDimensions[0]]
            columnNumber = plateColumnNumbers[col%plateDimensions[1]]
            key=-1
            if not colwrapBool:
                plateName = plateRowLetters[int(row/plateDimensions[0])]+str(int(col/plateDimensions[1])+1)
            else:
                plateName = plateRowLetters[0]+str(int((trueCol*colwrap*plateDimensions[1]+col)/plateDimensions[1])+1)
            if row not in rowLimits or col not in colLimits:            
                if row % plateDimensions[0] == 0 and col % plateDimensions[1] == 0:
                    infoList.append([rowLetter+str(columnNumber),plateName,plateName])
                else:
                    infoList.append([rowLetter+str(columnNumber),'DoNotLabel',plateName])
                pointList.append([x,y,key,int(row/plateDimensions[0]),int(col/plateDimensions[1]),0])
            else:
                infoList.append([rowLetter+str(columnNumber),'DoNotLabel','None'])
                pointList.append([x,y,key,int(row/plateDimensions[0]),int(col/plateDimensions[1]),-1])

    pointMatrix = np.matrix(pointList)
    infoMatrix = np.matrix(infoList)

    if colwrapBool:
        for i,matrix in enumerate([pointMatrix,infoMatrix]):
            newList = []
            for col in range(matrix.shape[1]):
                vals = np.reshape(matrix[:,col],(visualNumRowPlates*plateDimensions[0],visualNumColumnPlates*plateDimensions[1]))
                reshapedVals = blockshaped(np.asarray(vals), plateDimensions[0], plateDimensions[1])
                reshapedValList = [reshapedVals[x,:,:] for x in range(timepointPlateColumns)]
                #reshapedValList = [reshapedVals[x,:,:] for x in range(reshapedVals.shape[0])]
                reshapedMatrix = np.hstack(reshapedValList)
                newList.append(reshapedMatrix.flatten())
            if i == 0:
                pointMatrix = np.matrix(newList).T
            else:
                infoMatrix = np.matrix(newList).T
    
    pointDf = pd.DataFrame(pointMatrix,columns=['x','y','key','plateRow','plateColumn','blank'])
    infoDf = pd.DataFrame(infoMatrix,columns=['wellID','plateID','plateName'])
    
    #Make lines to separate plates
    vlineList = []
    for vline in range(1,visualNumColumnPlates):
        vlinex = (baseLayoutX[vline*plateDimensions[1]-1]+baseLayoutX[vline*plateDimensions[1]])/2.0
        vlineList.append(vlinex)
    hlineList = []
    for hline in range(1,visualNumRowPlates):
        hliney = (baseLayoutY[hline*plateDimensions[0]-1]+baseLayoutY[hline*plateDimensions[0]])/2.0
        hlineList.append(hliney)
    
    return pointDf,infoDf,vlineList,hlineList

class BlankSelectionPage(tk.Frame):
    def __init__(self, master,folderName,levels,levelValues,maxNumLevelValues,numPlates,plateDimensions,dtl,shp,bPage):
        numRowPlates = 1
        numColumnPlates = numPlates
        if numRowPlates == 1 and numColumnPlates > colwrap:
            colwrapBool = True
        else:
            colwrapBool = False

        self.root = master.root
        tk.Frame.__init__(self, master)
        
        global secondaryhomepage,dataTypeList,backPage
        secondaryhomepage = shp
        dataTypeList = dtl
        backPage = bPage

        tk.Label(self,text='Blank Selection Page:',font=('Arial',20)).grid(row=0,column=0)

        plotFrame = tk.Frame(self)
        plotFrame.grid(row=5,column=0)
        
        if plateDimensions[1] == 24:
            scalingFactor = 1.5
        else:
            scalingFactor = 1
        
        if colwrapBool:
            visualNumRowPlates = math.ceil(numColumnPlates/colwrap)
            visualNumColumnPlates = colwrap
        else:
            visualNumRowPlates = numRowPlates
            visualNumColumnPlates = numColumnPlates 
        fig = plt.figure(figsize=(3*visualNumColumnPlates*(plateDimensions[0]/8), 2*visualNumRowPlates*(plateDimensions[1]/12)*figLengthScaling),tight_layout=True)
        gs = fig.add_gridspec(1, 1)
        fig_ax1 = fig.add_subplot(gs[0])

        baseLayoutDf,infoDf,vlinelist,hlinelist = returnBaseLayout(plateDimensions,numRowPlates,numColumnPlates)

        self.currentLayout = baseLayoutDf.copy()
        def modifyPlatePlot():
            #Remove all axis elements
            fig_ax1.set_xlim((0, max(baseLayoutDf['x'])+1))
            fig_ax1.set_ylim((0, max(baseLayoutDf['y'])+1))
            fig_ax1.set_xticks([])
            fig_ax1.set_yticks([])
            fig_ax1.set_xlabel('')
            fig_ax1.set_ylabel('')
            #Add plate dividing lines
            for vlinex in vlinelist:
                fig_ax1.axvline(vlinex,linestyle=':',color='k')
            for hliney in hlinelist:
                fig_ax1.axhline(hliney,linestyle=':',color='k')
            #Add well IDs
            for row in range(baseLayoutDf.shape[0]):
                fig_ax1.annotate(infoDf.values[row][0],(baseLayoutDf.iloc[row,0],baseLayoutDf.iloc[row,1]),ha='center',va='center',size=7)
            #Add plateIDs
            for row in range(baseLayoutDf.shape[0]):
                if infoDf.values[row][1] != 'DoNotLabel':
                    fig_ax1.annotate(infoDf.values[row][1],(baseLayoutDf.iloc[row,0]-0.6,baseLayoutDf.iloc[row,1]+0.5),size=7,ha='center',va='center')
        
        self.canvas = FigureCanvasTkAgg(fig,master=plotFrame)
        modifyPlatePlot()
        self.path = sns.scatterplot(data=baseLayoutDf,x='x',y='y',ax=fig_ax1,color='#808080',s=200,marker='o',alpha=0.5)
        self.canvas.draw()
        self.background = self.canvas.copy_from_bbox(fig_ax1.bbox)
        self.trueBackground = self.canvas.copy_from_bbox(fig_ax1.bbox)

        #background = fig.canvas.copy_from_bbox(fig_ax1.bbox)
        #modifyPlatePlot()

        #self.canvas = FigureCanvasTkAgg(fig,master=plotFrame)
        #self.canvas.draw()
        self.canvas.get_tk_widget().pack()

        def line_select_callback(eclick, erelease):
            'eclick and erelease are the press and release events'
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata

        def toggle_selector(event):
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                toggle_selector.RS.set_active(True)
         
        rectpropsdict = {'facecolor':'#FF0000','alpha':0.2,'edgecolor':'#FF0000'}
        toggle_selector.RS = RectangleSelector(fig_ax1, line_select_callback,drawtype='box', useblit=True,button=[1, 3], minspanx=1, minspany=1,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
        self.ts = toggle_selector.RS
        self.canvas.mpl_connect('key_press_event', toggle_selector)
        
        def selectWells(mark):
            wellSelectionBox = toggle_selector.RS.corners
            ll = np.array([wellSelectionBox[0][0], wellSelectionBox[1][0]])  # lower-left
            ur = np.array([wellSelectionBox[0][2], wellSelectionBox[1][2]])  # upper-right
            inidx = np.all(np.logical_and(ll <= self.currentLayout.values[:,:2], self.currentLayout.values[:,:2] <= ur), axis=1)
            inbox = self.currentLayout.loc[inidx]
            changedInbox = inbox.copy()
            if mark:
                changedInbox['key'] = 0
                self.currentLayout.loc[inidx] = changedInbox 
                updatePlatePlot(changedInbox)
            else:
                changedInbox['key'] = -1
                self.currentLayout.loc[inidx] = changedInbox 
                updatePlatePlot(changedInbox)
        
        def updatePlatePlot(newlayout):
            self.path = sns.scatterplot(data=newlayout,x='x',y='y',ax=fig_ax1,s=200,marker='o',color=['#ffffff'])
            #g1 = sns.scatterplot(data=newlayout,x='x',y='y',ax=fig_ax1,s=200,style='key',color=['#ffffff'],markers=['o','X'],style_order=[-1,0])
            self.path = sns.scatterplot(data=newlayout,x='x',y='y',ax=fig_ax1,s=200,style='key',markers=['o','X'],style_order=[-1,0],alpha=0.5,hue_order=[-1,0],hue='key',palette=['#808080','#FF0000'])
            fig_ax1.legend_.remove()
            fig_ax1.draw_artist(self.path)
            self.canvas.blit(fig_ax1.bbox)
            self.background = self.canvas.copy_from_bbox(fig_ax1.bbox)
            toggle_selector.RS.background = self.background 
            #self.canvas.draw()
        
        def clearWells():
            MsgBox = tk.messagebox.askquestion ('Warning','Are you sure you want to clear all wells?',icon = 'warning')
            if MsgBox == 'yes':
                self.currentLayout = baseLayoutDf.copy()
                self.path = sns.scatterplot(data=self.currentLayout,x='x',y='y',ax=fig_ax1,s=200,style='key',color=['#ffffff'],markers=['o','X'],style_order=[-1,0])
                self.path = sns.scatterplot(data=self.currentLayout,x='x',y='y',ax=fig_ax1,s=200,marker='o',alpha=0.5,color=['#808080'])
                fig_ax1.legend_.remove()
                self.canvas.restore_region(self.trueBackground)
                self.canvas.blit(fig_ax1.bbox)
                self.background = self.canvas.copy_from_bbox(fig_ax1.bbox)
                toggle_selector.RS.background = self.background 
                #self.currentLayout = baseLayoutDf.copy()
                #self.canvas.draw()
        
        def collectInputs():
            master.switch_frame(PlateLayoutPage,folderName,self.currentLayout['key'],levels,levelValues,maxNumLevelValues,numRowPlates,numColumnPlates,plateDimensions,dataTypeList)

        def updateExperimentPlot():
            self.path = sns.scatterplot(data=self.currentLayout,x='x',y='y',ax=fig_ax1,s=200,marker='o',color=['#ffffff'])
            #g1 = sns.scatterplot(data=newlayout,x='x',y='y',ax=fig_ax1,s=200,style='key',color=['#ffffff'],markers=['o','X'],style_order=[-1,0])
            self.path = sns.scatterplot(data=self.currentLayout,x='x',y='y',ax=fig_ax1,s=200,style='key',markers=['o','X'],style_order=[-1,0],alpha=0.5,hue_order=[-1,0],hue='key',palette=['#808080','#FF0000'])
            fig_ax1.legend_.remove()
            fig_ax1.draw_artist(self.path)
            self.canvas.blit(fig_ax1.bbox)
            self.background = self.canvas.copy_from_bbox(fig_ax1.bbox)
            toggle_selector.RS.background = self.background 
            #self.canvas.draw()

        def repeatSelection(direction):
            reshapedList = []
            for column in self.currentLayout.columns:
                reshapedList.append(np.matrix(np.reshape(self.currentLayout[column].values,(plateDimensions[0]*numRowPlates,plateDimensions[1]*numColumnPlates))))
            keyMatrix = reshapedList[2]
            
            currentPlateElementsRow = np.where(keyMatrix != -1)[0]
            currentPlateElementsCol = np.where(keyMatrix != -1)[1]
            plateToRepeat = keyMatrix[np.ix_(np.unique(currentPlateElementsRow),np.unique(currentPlateElementsCol))]
            fullKeyLayout = keyMatrix.copy()
            if direction == 'h':
                repeatedLayout = np.tile(plateToRepeat,numColumnPlates)
                fullKeyLayout[np.ix_(np.unique(currentPlateElementsRow),list(range(plateDimensions[1]*numColumnPlates)))] = repeatedLayout
            else:
                repeatedLayout = np.tile(plateToRepeat,(numRowPlates,1))
                fullKeyLayout[np.ix_(list(range(plateDimensions[0]*numRowPlates)),np.unique(currentPlateElementsCol))] = repeatedLayout
            
            unrolledKeys = fullKeyLayout.flatten()
            self.currentLayout['key'] = unrolledKeys.T
            updateExperimentPlot()

        selectionWindow = tk.Frame(self)
        selectionWindow.grid(row=1,column=0)
        tk.Button(selectionWindow, text="Deselect Wells",command=lambda: selectWells(False)).grid(row=0,column=0)
        tk.Button(selectionWindow, text="Select Wells",command=lambda: selectWells(True)).grid(row=0,column=1)
        
        
        actionWindow = tk.Frame(self)
        actionWindow.grid(row=2,column=0)
        tk.Button(actionWindow, text="Repeat Plate Horizontally",command=lambda: repeatSelection('h')).grid(row=0,column=0)
        tk.Button(actionWindow, text="Repeat Plate Vertically",command=lambda: repeatSelection('v')).grid(row=0,column=1)

        tk.Button(self, text="Clear All Wells",command=lambda: clearWells()).grid(row=3,column=0)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.grid(row=4,column=0)
        
        tk.Button(buttonWindow, text="OK",command=lambda: collectInputs(),font='Helvetica 14 bold').grid(row=0,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(secondaryhomepage,folderName,bPage)).grid(row=0,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=0,column=2)

def createLayoutVisual(baseLayoutDf,currentLayout,levelIndex,currentLevel,levelValues,plateDimensions,numRowPlates,numColumnPlates,dt,infoDf,vlinelist,hlinelist):
    
    maxTextLen = len(str(currentLevel))
    for levelVal in levelValues[levelIndex]:
        if len(str(levelVal)) > maxTextLen:
            maxTextLen = len(str(levelVal))
    
    if numRowPlates == 1 and numColumnPlates > colwrap:
        colwrapBool = True
    else:
        colwrapBool = False

    if colwrapBool:
        visualNumRowPlates = math.ceil(numColumnPlates/colwrap)
        visualNumColumnPlates = colwrap
    else:
        visualNumRowPlates = numRowPlates
        visualNumColumnPlates = numColumnPlates 
    fig = plt.figure(figsize=(3*visualNumColumnPlates*(plateDimensions[0]/8)+1+(0.125*maxTextLen), 2.5*visualNumRowPlates*(plateDimensions[1]/12)*figLengthScaling),tight_layout=True)
    #fig_ax1 = fig.add_subplot(111)
    gs = fig.add_gridspec(1, 1)
    fig_ax1 = fig.add_subplot(gs[0])
    
    #Remove all axis elements
    fig_ax1.set_xlim((0, max(baseLayoutDf['x'])+1))
    fig_ax1.set_ylim((0, max(baseLayoutDf['y'])+1))
    fig_ax1.set_xticks([])
    fig_ax1.set_yticks([])
    fig_ax1.set_xlabel('')
    fig_ax1.set_ylabel('')
    #Add plate dividing lines
    for vlinex in vlinelist:
        fig_ax1.axvline(vlinex,linestyle=':',color='k')
    for hliney in hlinelist:
        fig_ax1.axhline(hliney,linestyle=':',color='k')
    #Add well IDs
    for row in range(baseLayoutDf.shape[0]):
        fig_ax1.annotate(infoDf.values[row][0],(baseLayoutDf.iloc[row,0],baseLayoutDf.iloc[row,1]),ha='center',va='center',size=7)
    #Add plateIDs
    for row in range(baseLayoutDf.shape[0]):
        if infoDf.values[row][1] != 'DoNotLabel':
            fig_ax1.annotate(infoDf.values[row][1],(baseLayoutDf.iloc[row,0]-0.6,baseLayoutDf.iloc[row,1]+0.5),size=7,ha='center',va='center')
    
    currentpalette = ['#808080']+sns.color_palette('husl',len(levelValues[levelIndex])).as_hex()
    newLayoutDf = baseLayoutDf.copy()
    newLayoutDf['key'] = currentLayout.flatten()

    g1 = sns.scatterplot(data=baseLayoutDf,x='x',y='y',ax=fig_ax1,color='#ffffff',s=200,marker='o')
    if -1 in list(newLayoutDf['key']) or 'blank' in list(newLayoutDf['key']):
        hueorder = [-1]
        modifiedPalette = ['#808080']
        ogpalette = modifiedPalette.copy()
    else:
        hueorder = []
        modifiedPalette = []
        ogpalette = modifiedPalette.copy()
    for i in range(len(levelValues[levelIndex])):
        if i in list(pd.unique(newLayoutDf['key'])):
            hueorder.append(i)
            modifiedPalette.append(currentpalette[i+1])
    blanksExist = False 
    for i,val in enumerate(newLayoutDf['key']):
        if val == 'blank':
            newLayoutDf['key'].iloc[i] = -1
            newLayoutDf['blank'].iloc[i] = -1
            blanksExist = True
        else:
            newLayoutDf['blank'].iloc[i] = 0 
    g1 = sns.scatterplot(data=newLayoutDf,x='x',y='y',ax=fig_ax1,hue='key',hue_order=hueorder,palette=modifiedPalette,s=200,markers=['X','o'],alpha=0.5,style='blank',style_order=[-1,0])
    titlelabels = [currentLevel+': ']
    titlecolors = ['black']
    
    for i,levelValue in enumerate(levelValues[levelIndex]):
        titlelabels.append(str(levelValue))
        titlecolors.append(modifiedPalette[i+len(ogpalette)])
    #plt.savefig('plateLayout-'+currentLevel+'.png',bbox_inches='tight')
    if 'plateLayouts' not in os.listdir('plots'):
        subprocess.run(['mkdir','plots/plateLayouts'])
    legendHandlesLabels = fig_ax1.get_legend_handles_labels()
    i=0
    newLegendLabels = []
    newLegendHandles = []
    legendHandles = legendHandlesLabels[0]
    currentLevelValues = levelValues[levelIndex]
    if blanksExist:
        offset = 1
    else:
        offset = 0
    for legendHandle in legendHandles:
        #Skips the style legend handles
        if i < len(currentLevelValues)+1:
            if i == 0:
                #modifiedLevel = currentLevel.translate(str.maketrans({"-":  r"\-","]":  r"\]","\\": r"\\","^":  r"\^","$":  r"\$","*":  r"\*",".":  r"\.","_":  r"\_","(":  r"\(",")":  r"\)","[":  r"\[","%": r"\%"}))
                modifiedLevel = currentLevel
                newLegendLabels.append('$\\bf'+modifiedLevel+'$')
                newLegendHandles.append(legendHandle)
            else:
                newLegendLabels.append(currentLevelValues[i-1])
                newLegendHandles.append(legendHandles[i+offset])
        i+=1
    fig_ax1.legend(bbox_to_anchor=(1, 1),frameon=False,handles=newLegendHandles, labels=newLegendLabels)
    plt.savefig('plots/plateLayouts/plateLayout-'+currentLevel+'-'+dt+'.png',bbox_inches='tight')
    plt.clf()

class PlateLayoutPage(tk.Frame):
    def __init__(self, master,folderName,blankWells,levels,levelValues,maxNumLevelValues,numRowPlates,numColumnPlates,plateDimensions,dataTypeList):
        
        self.root = master.root
        tk.Frame.__init__(self, master)
        
        tk.Label(self,text='Label Layout Page:',font=('Arial',20)).grid(row=0,column=0)
        
        plotFrame = tk.Frame(self)
        plotFrame.grid(row=9,column=0)
        
        if numRowPlates == 1 and numColumnPlates > colwrap:
            colwrapBool = True
        else:
            colwrapBool = False
        
        if colwrapBool:
            visualNumRowPlates = math.ceil(numColumnPlates/colwrap)
            visualNumColumnPlates = colwrap
        else:
            visualNumRowPlates = numRowPlates
            visualNumColumnPlates = numColumnPlates 
        fig = plt.figure(figsize=(3*visualNumColumnPlates*(plateDimensions[0]/8), 2*visualNumRowPlates*(plateDimensions[1]/12)*figLengthScaling),tight_layout=True)
        #fig = plt.figure(figsize=(3*numColumnPlates*(plateDimensions[0]/8), 2*numRowPlates*(plateDimensions[1]/12)*figLengthScaling),tight_layout=True)
        gs = fig.add_gridspec(1, 1)
        fig_ax1 = fig.add_subplot(gs[0])
        

        #baseLayoutDf,infoDf,vlinelist,hlinelist = returnBaseLayout([16,24],2,2)
        baseLayoutDf,infoDf,vlinelist,hlinelist = returnBaseLayout(plateDimensions,numRowPlates,numColumnPlates)
        baseLayoutDf['blank'] = blankWells
        
        blankMatrix = np.matrix(np.reshape(baseLayoutDf['blank'].values,(plateDimensions[0]*numRowPlates,plateDimensions[1]*numColumnPlates)))
        self.blankRowsToDisregard = []
        self.blankColumnsToDisregard = []
        self.blankPlateRowsToDisregard = []
        self.blankPlateColumnsToDisregard = []
        for row in range(blankMatrix.shape[0]):
            if np.all(blankMatrix[row,:] == 0):
                self.blankRowsToDisregard.append(row)
                if row < plateDimensions[0]:
                    self.blankPlateRowsToDisregard.append(row)
        for col in range(blankMatrix.shape[1]):
            if np.all(blankMatrix[:,col] == 0):
                self.blankColumnsToDisregard.append(col)
                if col < plateDimensions[1]:
                    self.blankPlateColumnsToDisregard.append(col)
        
        self.currentLayout = baseLayoutDf.copy()
        self.infoDf = infoDf.copy()
        self.allLayouts = [baseLayoutDf.copy()]*len(levels)
        self.levelIndex = 0
        self.levelValueIndex = 0
        self.currentpalette = ['#808080']+sns.color_palette('husl',len(levelValues[self.levelIndex])).as_hex()

        def modifyPlatePlot():
            #Remove all axis elements
            fig_ax1.set_xlim((0, max(baseLayoutDf['x'])+1))
            fig_ax1.set_ylim((0, max(baseLayoutDf['y'])+1))
            fig_ax1.set_xticks([])
            fig_ax1.set_yticks([])
            fig_ax1.set_xlabel('')
            fig_ax1.set_ylabel('')
            #Add plate dividing lines
            for vlinex in vlinelist:
                fig_ax1.axvline(vlinex,linestyle=':',color='k')
            for hliney in hlinelist:
                fig_ax1.axhline(hliney,linestyle=':',color='k')
            #Add well IDs
            for row in range(baseLayoutDf.shape[0]):
                fig_ax1.annotate(infoDf.values[row][0],(baseLayoutDf.iloc[row,0],baseLayoutDf.iloc[row,1]),ha='center',va='center',size=7)
            #Add plateIDs
            for row in range(baseLayoutDf.shape[0]):
                if infoDf.values[row][1] != 'DoNotLabel':
                    fig_ax1.annotate(infoDf.values[row][1],(baseLayoutDf.iloc[row,0]-0.6,baseLayoutDf.iloc[row,1]+0.5),size=7,ha='center',va='center')

        self.canvas = FigureCanvasTkAgg(fig,master=plotFrame)
        modifyPlatePlot()
        self.path = sns.scatterplot(data=baseLayoutDf,x='x',y='y',ax=fig_ax1,color='#808080',alpha=0.5,s=200,markers=['o','X'],style='blank',style_order=[-1,0])
        fig_ax1.legend_.remove()
        self.canvas.draw()
        self.background = self.canvas.copy_from_bbox(fig_ax1.bbox)
        self.trueBackground = self.canvas.copy_from_bbox(fig_ax1.bbox)
        
        #self.canvas.restore_region(background)
        #fig_ax1.draw_artist(self.path)
        #self.canvas.blit(fig_ax1.bbox)

        self.canvas.get_tk_widget().pack()
        
        def line_select_callback(eclick, erelease):
            'eclick and erelease are the press and release events'
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata

        def toggle_selector(event):
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                toggle_selector.RS.set_active(True)
         
        rectpropsdict = {'facecolor':self.currentpalette[self.levelValueIndex+1],'alpha':0.2,'edgecolor':self.currentpalette[self.levelValueIndex+1]}
        toggle_selector.RS = RectangleSelector(fig_ax1, line_select_callback,drawtype='box', useblit=True,button=[1, 3], minspanx=1, minspany=1,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
        self.ts = toggle_selector.RS
        toggle_selector.RS.background = self.trueBackground 

        def updatePlatePlot(newlayout,key):
            currentcolor = self.currentpalette[key]
            self.path = sns.scatterplot(data=newlayout,x='x',y='y',ax=fig_ax1,s=200,color='#ffffff',marker='o')
            self.path = sns.scatterplot(data=newlayout,x='x',y='y',ax=fig_ax1,s=200,alpha=0.5,markers=['o','X'],style='blank',style_order=[-1,0],color=currentcolor)
            fig_ax1.legend_.remove()
        
        def clearWells():
            MsgBox = tk.messagebox.askquestion ('Warning','Are you sure you want to clear all wells?',icon = 'warning')
            if MsgBox == 'yes':
                self.currentLayout = baseLayoutDf.copy()
                self.path = sns.scatterplot(data=self.currentLayout,x='x',y='y',ax=fig_ax1,s=200,markers=['o','X'],style='blank',style_order=[-1,0],color='#ffffff')
                self.path = sns.scatterplot(data=self.currentLayout,x='x',y='y',ax=fig_ax1,s=200,markers=['o','X'],alpha=0.5,style='blank',style_order=[-1,0],color='#808080')
                fig_ax1.legend_.remove()
                self.canvas.restore_region(self.trueBackground)
                self.canvas.blit(fig_ax1.bbox)
                self.background = self.trueBackground
                toggle_selector.RS.background = self.trueBackground 
                #self.canvas.draw()

        def changeLevelValue(advance):
            if advance:
                self.levelValueIndex += 1
            else:
                self.levelValueIndex -= 1
            
            if self.levelValueIndex == len(levelValues[self.levelIndex])-1:
                self.nextLevelValueButton['state'] = 'disabled'
            else:
                self.nextLevelValueButton['state'] = 'normal'
            if self.levelValueIndex == 0:
                self.previousLevelValueButton['state'] = 'disabled'
            else:
                self.previousLevelValueButton['state'] = 'normal'
            
            self.levelValueIndex = max([0,self.levelValueIndex])
            self.levelValueIndex = min([len(levelValues[self.levelIndex])-1,self.levelValueIndex])
            
            changeLevelValueInLevelValueLabelList()
            rectpropsdict = {'facecolor':self.currentpalette[self.levelValueIndex+1],'alpha':0.2,'edgecolor':self.currentpalette[self.levelValueIndex+1]}
            toggle_selector.RS = RectangleSelector(fig_ax1, line_select_callback, useblit=True,drawtype='box',button=[1, 3], minspanx=1, minspany=1,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
            fig_ax1.draw_artist(self.path)
            self.canvas.blit(fig_ax1.bbox)
            self.background = self.canvas.copy_from_bbox(fig_ax1.bbox)
            toggle_selector.RS.background = self.background 
            #self.canvas.draw()
        
        def updateExperimentPlot():
            self.path = sns.scatterplot(data=baseLayoutDf,x='x',y='y',ax=fig_ax1,color='#ffffff',s=200,marker='o')
            if -1 in list(self.currentLayout['key']):
                hueorder = [-1]
                modifiedPalette = ['#808080']
            else:
                hueorder = []
                modifiedPalette = []
            for i in range(len(levelValues[self.levelIndex])):
                if i in list(pd.unique(self.currentLayout['key'])):
                    hueorder.append(i)
                    modifiedPalette.append(self.currentpalette[i+1])
            self.path = sns.scatterplot(data=self.currentLayout,x='x',y='y',ax=fig_ax1,hue='key',hue_order=hueorder,palette=modifiedPalette,s=200,markers=['o','X'],alpha=0.5,style='blank',style_order=[-1,0])
            fig_ax1.legend_.remove()
            toggle_selector.RS = RectangleSelector(fig_ax1, line_select_callback, useblit=True,drawtype='box',button=[1, 3], minspanx=1, minspany=1,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
            fig_ax1.draw_artist(self.path)
            self.canvas.blit(fig_ax1.bbox)
            self.background = self.canvas.copy_from_bbox(fig_ax1.bbox)
            toggle_selector.RS.background = self.background 
            #self.canvas.draw()
        
        def changeLevel(advance):
            self.levelValueIndex = 0
            self.allLayouts[self.levelIndex] = self.currentLayout.copy()
            if advance:
                self.levelIndex+=1
            else:
                self.levelIndex-=1
            
            if self.levelIndex == len(levels)-1:
                self.FinishButton['state'] = 'normal'
                self.nextLevelButton['state'] = 'disabled'
            else:
                self.FinishButton['state'] = 'disabled'
                self.nextLevelButton['state'] = 'normal'
            if self.levelIndex == 0:
                self.previousLevelButton['state'] = 'disabled'
            else:
                self.previousLevelButton['state'] = 'normal'
            
            self.nextLevelValueButton['state'] = 'normal'
            self.previousLevelValueButton['state'] = 'disabled'
            
            self.levelIndex = max([0,self.levelIndex])
            self.levelIndex = min([len(levels)-1,self.levelIndex])
            
            self.currentLayout = self.allLayouts[self.levelIndex].copy()
            self.currentpalette = ['#808080']+sns.color_palette('husl',len(levelValues[self.levelIndex])).as_hex()
            
            changeLevelInLevelLabelList()
            changeLevelInLevelValueLabelList()
            #rectpropsdict = {'facecolor':self.currentpalette[self.levelValueIndex+1],'alpha':0.2,'edgecolor':self.currentpalette[self.levelValueIndex+1]}
            #toggle_selector.RS = RectangleSelector(fig_ax1, line_select_callback,drawtype='box', useblit=True,button=[1, 3], minspanx=1, minspany=1,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
            updateExperimentPlot()

        #Level labels
        levelLabelWindow = tk.Frame(self)
        levelLabelWindow.grid(row=1,column=0)
        
        levelTitle = tk.Label(levelLabelWindow,text='Level: ')
        levelTitle.grid(row=0,column=0,sticky=tk.W)
        levelTitle.configure(font='Helvetica 18 bold')
        self.levelLabelList = []
        for i,level in enumerate(levels):
            currentlabel = tk.Label(levelLabelWindow,text=level)
            if i == self.levelIndex:
                currentlabel.configure(font='Helvetica 18',borderwidth=2, relief="solid")
            else:
                currentlabel.configure(font='Helvetica 18',borderwidth=0, relief="solid")
            currentlabel.grid(row=0,column=i+1,sticky=tk.W, padx=2)
            self.levelLabelList.append(currentlabel)
        
        def changeLevelInLevelLabelList():
            for i,level in enumerate(levels):
                currentlabel = self.levelLabelList[i]
                if i == self.levelIndex:
                    currentlabel.configure(font='Helvetica 18',borderwidth=2, relief="solid")
                else:
                    currentlabel.configure(font='Helvetica 18',borderwidth=0, relief="solid")
        
        self.previousLevelButton = tk.Button(levelLabelWindow, text="Previous",command=lambda: changeLevel(False))
        self.previousLevelButton.grid(row=0,column=len(levels)+1,sticky=tk.W)
        self.previousLevelButton['state'] = 'disabled'
        self.nextLevelButton = tk.Button(levelLabelWindow, text="Next",command=lambda: changeLevel(True))
        self.nextLevelButton.grid(row=0,column=len(levels)+2,sticky=tk.W)

        #Level value labels
        levelValueLabelWindow = tk.Frame(self)
        levelValueLabelWindow.grid(row=2,column=0)

        levelValueTitle = tk.Label(levelValueLabelWindow,text='Level Value: ')
        levelValueTitle.grid(row=0,column=0,sticky=tk.W+tk.E)
        levelValueTitle.configure(font='Helvetica 18 bold')
        #Need to create maxNumLevelValue labels; only fill in appropriate labels for current level
        self.levelValueLabelList = []
        maxNumLevelValues = len(max(levelValues,key=len))
        if maxNumLevelValues%2 == 0:
            maxNumLevelValues+=1
        labelIndexList = list(range(maxNumLevelValues))
        median = labelIndexList[int(len(labelIndexList)/2)]
        
        self.blanktext = ''

        #Instantiate all level value labels
        numCurrentLevelValues = len(levelValues[self.levelIndex])
        currentLabelIndexList = labelIndexList[median-int(numCurrentLevelValues/2):median-int(numCurrentLevelValues/2)+numCurrentLevelValues]
        j=0
        for i in range(maxNumLevelValues):
            if i in currentLabelIndexList:
                levelValue = levelValues[self.levelIndex][j]
                currentlabel = tk.Label(levelValueLabelWindow,text=levelValue,fg = self.currentpalette[j+1])
                if j == self.levelValueIndex:
                    currentlabel.configure(font='Helvetica 18',borderwidth=2, relief="solid")
                else:
                    currentlabel.configure(font='Helvetica 18',borderwidth=0, relief="solid")
                j+=1
            else:
                currentlabel = tk.Label(levelValueLabelWindow,text=self.blanktext)
            currentlabel.grid(row=0,column=i+1,sticky=tk.W+tk.E,padx=2)
            self.levelValueLabelList.append(currentlabel)
        
        def changeLevelValueInLevelValueLabelList():
            numCurrentLevelValues = len(levelValues[self.levelIndex])
            currentLabelIndexList = labelIndexList[median-int(numCurrentLevelValues/2):median-int(numCurrentLevelValues/2)+numCurrentLevelValues]
            j=0
            for i in range(maxNumLevelValues):
                currentlabel = self.levelValueLabelList[i]
                if i in currentLabelIndexList:
                    if currentlabel['text'] != self.blanktext:
                        if j == self.levelValueIndex:
                            currentlabel.configure(font='Helvetica 18',borderwidth=2, relief="solid")
                        else:
                            currentlabel.configure(font='Helvetica 18',borderwidth=0, relief="solid")
                    else:
                        currentlabel.configure(font='Helvetica 18',borderwidth=0, relief="solid")
                    j+=1

        self.previousLevelValueButton = tk.Button(levelValueLabelWindow, text="Previous",command=lambda: changeLevelValue(False))
        self.previousLevelValueButton.grid(row=0,column=maxNumLevelValues+1,sticky=tk.W)
        self.previousLevelValueButton['state'] = 'disabled'
        self.nextLevelValueButton = tk.Button(levelValueLabelWindow, text="Next",command=lambda: changeLevelValue(True))
        self.nextLevelValueButton.grid(row=0,column=maxNumLevelValues+2,sticky=tk.W)
        
        def changeLevelInLevelValueLabelList():
            numCurrentLevelValues = len(levelValues[self.levelIndex])
            currentLabelIndexList = labelIndexList[median-int(numCurrentLevelValues/2):median-int(numCurrentLevelValues/2)+numCurrentLevelValues]
            j = 0
            for i in range(maxNumLevelValues):
                currentlabel = self.levelValueLabelList[i]
                if i in currentLabelIndexList:
                    levelValue = levelValues[self.levelIndex][j]
                    currentlabel.configure(text=levelValue,fg = self.currentpalette[j+1])
                    if j == self.levelValueIndex:
                        currentlabel.configure(font='Helvetica 18',borderwidth=2, relief="solid")
                    else:
                        currentlabel.configure(font='Helvetica 18',borderwidth=0, relief="solid")
                    j+=1
                else:
                    currentlabel.configure(text=self.blanktext,borderwidth=0, relief="solid")
            if numColumnPlates <= 2:
                currentlabel.update()
        
        def selectWells(mark,immediate):

            wellSelectionBox = toggle_selector.RS.corners
            ll = np.array([wellSelectionBox[0][0], wellSelectionBox[1][0]])  # lower-left
            ur = np.array([wellSelectionBox[0][2], wellSelectionBox[1][2]])  # upper-right
            inidx = np.all(np.logical_and(ll <= self.currentLayout.values[:,:2], self.currentLayout.values[:,:2] <= ur), axis=1)
            inbox = self.currentLayout.loc[inidx]
            changedInbox = inbox.copy().loc[inbox.blank != 0,:]
            #Remove blank wells from selection to change
            inidx2 = inidx.copy()
            for row in range(self.currentLayout.shape[0]):
                if self.currentLayout.iloc[row,:]['blank'] == 0:
                    inidx2[row] = False

            if mark:
                changedInbox.loc[:,'key'] = self.levelValueIndex
                self.currentLayout.loc[inidx2] = changedInbox
                updatePlatePlot(changedInbox,self.levelValueIndex+1)
                if immediate == 'yes' and self.levelValueIndex != maxNumLevelValues-1:
                    changeLevelValue(True) 
                else:
                    rectpropsdict = {'facecolor':self.currentpalette[self.levelValueIndex+1],'alpha':0.2,'edgecolor':self.currentpalette[self.levelValueIndex+1]}
                    toggle_selector.RS = RectangleSelector(fig_ax1, line_select_callback, useblit=True,drawtype='box',button=[1, 3], minspanx=1, minspany=1,spancoords='pixels',interactive=True,rectprops=rectpropsdict)
                    fig_ax1.draw_artist(self.path)
                    self.canvas.blit(fig_ax1.bbox)
                    self.background = self.canvas.copy_from_bbox(fig_ax1.bbox)
                    toggle_selector.RS.background = self.background 
            else:
                changedInbox['key'] = -1
                self.currentLayout.loc[inidx] = changedInbox 
                updatePlatePlot(changedInbox,0)
                toggle_selector.RS.background = self.background 
                fig_ax1.draw_artist(self.path)
                self.canvas.blit(fig_ax1.bbox)

        def copySelection():
            wellSelectionBox = toggle_selector.RS.corners
            ll = np.array([wellSelectionBox[0][0], wellSelectionBox[1][0]])  # lower-left
            ur = np.array([wellSelectionBox[0][2], wellSelectionBox[1][2]])  # upper-right
            inidx = np.all(np.logical_and(ll <= self.currentLayout.values[:,:2], self.currentLayout.values[:,:2] <= ur), axis=1)
            inbox = self.currentLayout.loc[inidx]
            changedInbox = inbox.copy().loc[inbox.blank != 0,:]
            #Remove blank wells from selection to change
            self.storedSelection = changedInbox['key'].values
        
        def pasteSelection():
            wellSelectionBox = toggle_selector.RS.corners
            ll = np.array([wellSelectionBox[0][0], wellSelectionBox[1][0]])  # lower-left
            ur = np.array([wellSelectionBox[0][2], wellSelectionBox[1][2]])  # upper-right
            inidx = np.all(np.logical_and(ll <= self.currentLayout.values[:,:2], self.currentLayout.values[:,:2] <= ur), axis=1)
            inbox = self.currentLayout.loc[inidx]
            changedInbox = inbox.copy().loc[inbox.blank != 0,:]
            #Remove blank wells from selection to change
            inidx2 = inidx.copy()
            for row in range(self.currentLayout.shape[0]):
                if self.currentLayout.iloc[row,:]['blank'] == 0:
                    inidx2[row] = False
            changedInbox.loc[:,'key'] = self.storedSelection
            self.currentLayout.loc[inidx2] = changedInbox
            updateExperimentPlot()

        def operateOnSelection(action,area,direction):
            
            keyMatrix = np.matrix(np.reshape(self.currentLayout['key'].values,(plateDimensions[0]*numRowPlates,plateDimensions[1]*numColumnPlates)))
            selectedElements = np.where(keyMatrix != -1)
            selectionToOperateOn = keyMatrix[np.ix_(np.unique(selectedElements[0]),np.unique(selectedElements[1]))]
            offset = [1,1]
            if area == 'experiment':
                areashape = keyMatrix.shape
                blankRows = self.blankRowsToDisregard
                blankColumns = self.blankColumnsToDisregard
            else:
                areashape = (plateDimensions[0],plateDimensions[1])
                blankRows = self.blankPlateRowsToDisregard
                blankColumns = self.blankPlateColumnsToDisregard
            blankArray = [blankRows,blankColumns]
            tiledLevelLayout = np.ones(keyMatrix.shape)*-1
            tilingMatrixList = []
            if action == 'tile':
                if direction == 'h':
                    numTiles = math.ceil(len(levelValues[self.levelIndex])/(areashape[1]-len(blankColumns)))
                    spillover = (areashape[1]-len(blankColumns))%len(levelValues[self.levelIndex])
                    shapeIndex = 1 
                else:
                    numTiles = math.ceil(len(levelValues[self.levelIndex])/(areashape[0]-len(blankRows)))
                    spillover = (areashape[0]-len(blankRows))%len(levelValues[self.levelIndex])
                    shapeIndex = 0
                if spillover == 0:
                    for tileIndex in range(len(levelValues[self.levelIndex])):
                        tilingMatrix = np.ones(selectionToOperateOn.shape) * tileIndex
                        tilingMatrixList.append(tilingMatrix)
                    if direction == 'h':
                        fullTilingMatrix = np.hstack(tilingMatrixList)
                    else:
                        fullTilingMatrix = np.vstack(tilingMatrixList)
                else:
                    splitLevelValues = np.array_split(np.array(levelValues[self.levelIndex]),int(numTiles))
                    splitLevelValues2 = np.array_split(levelValues[self.levelIndex],numTiles+1)
                    fullTilingMatrixList = []
                    for i,splitLevelValue in enumerate(splitLevelValues):
                        partialTilingMatrixList = []
                        for tileIndex in range(len(splitLevelValue)):
                            tilingMatrix = np.ones(selectionToOperateOn.shape) * (tileIndex + i*len(splitLevelValues[0]))
                            partialTilingMatrixList.append(tilingMatrix)
                        if len(levelValues[self.levelIndex])%(areashape[shapeIndex]-len(blankArray[shapeIndex])) != 0:
                            for spilloverValue in range(int(spillover/selectionToOperateOn.shape[shapeIndex])):
                                tilingMatrix = np.ones(selectionToOperateOn.shape) * -1
                                partialTilingMatrixList.append(tilingMatrix)
                        if direction == 'h':
                            partialTilingMatrix = np.hstack(partialTilingMatrixList)
                        else:
                            partialTilingMatrix = np.vstack(partialTilingMatrixList)
                        fullTilingMatrixList.append(partialTilingMatrix)
                    fullTilingMatrix = np.ones(areashape)*-1
                    if direction == 'h':
                        fullTilingMatrix[:int(numTiles),:int((areashape[1]-len(blankColumns)))] = np.vstack(fullTilingMatrixList)
                    else:
                        fullTilingMatrix[:int((areashape[0]-len(blankRows))),:int(numTiles)] = np.hstack(fullTilingMatrixList)
            else:
                if direction == 'h':
                    numRepeats = int((areashape[1]-len(self.blankColumnsToDisregard))/selectionToOperateOn.shape[1])
                    fullTilingMatrix = np.tile(selectionToOperateOn,numRepeats)
                else:
                    numRepeats = int((areashape[0]-len(self.blankRowsToDisregard))/selectionToOperateOn.shape[0])
                    fullTilingMatrix = np.tile(selectionToOperateOn,(numRepeats,1))
            tilerow = 0
            if area == 'plate':
                if direction == 'h':
                    shapeIndex = 0
                else:
                    shapeIndex = 1
                #Not sure what this does
                if fullTilingMatrix.shape[shapeIndex] != (areashape[shapeIndex]-len(blankArray[shapeIndex])):
                    off = math.ceil((areashape[shapeIndex]-len(blankArray[shapeIndex]))/fullTilingMatrix.shape[shapeIndex])
                    offset[shapeIndex] = off
            if action == 'repeat':
                for row in range(int(areashape[0]/offset[0])):
                    tilecol = 0
                    for col in range(int(areashape[1]/offset[1])):
                        if row not in blankArray[0]:
                            if col not in blankArray[1]:
                                tiledLevelLayout[row,col] = fullTilingMatrix[tilerow,tilecol]
                                tilecol += 1
                    if row not in blankArray[0]:
                        tilerow += 1
            else:
                for row in range(areashape[0]):
                    tilecol = 0
                    for col in range(areashape[1]):
                        if row not in blankArray[0]:
                            if col not in blankArray[1]:
                                tiledLevelLayout[row,col] = fullTilingMatrix[tilerow,tilecol]
                                tilecol += 1
                                if tilecol >= fullTilingMatrix.shape[1]:
                                    break
                    if row not in blankArray[0]:
                        tilerow += 1
                        if tilerow >= fullTilingMatrix.shape[0]:
                            break
            fullTiledLevelLayout = tiledLevelLayout.copy()
            """
            for row in range(int(areashape[0]/offset[0])):
                tilecol = 0
                for col in range(int(areashape[1]/offset[1])):
                    if row not in blankArray[0]:
                        print('wat2')
                        if col not in blankArray[1]:
                            print('wat3')
                            tiledLevelLayout[row,col] = fullTilingMatrix[tilerow,tilecol]
                            tilecol += 1
                if row not in blankArray[0]:
                    tilerow += 1
            fullTiledLevelLayout = tiledLevelLayout.copy()
            """
            unrolledKeys = fullTiledLevelLayout.flatten()
            if action == 'tile':
                self.currentLayout['key'] = unrolledKeys
            else:
                self.currentLayout['key'] = unrolledKeys.T
            updateExperimentPlot()

        def collectInputs():
            if self.levelIndex == len(levels)-1:
                self.allLayouts[self.levelIndex] = self.currentLayout.copy()
                finalLayoutDict = {}
                levelValueDict = {}
                keyDict = {}
                for levelKey,layout in enumerate(self.allLayouts):
                    keyMatrix = np.reshape(layout['key'].values,(plateDimensions[0]*numRowPlates,plateDimensions[1]*numColumnPlates))
                    levelValueMatrix = np.empty(keyMatrix.shape, dtype=object)
                    for row in range(keyMatrix.shape[0]):
                        for col in range(keyMatrix.shape[1]):
                            levelValueKey = int(keyMatrix[row][col])
                            if levelValueKey == -1:
                                levelValueMatrix[row,col] = 'blank'
                            else:
                                levelValueMatrix[row,col] = levelValueKey

                    levelValueDict[levelKey] = levelValueMatrix

                wellIDMatrix = np.reshape(self.infoDf['wellID'].values,(plateDimensions[0]*numRowPlates,plateDimensions[1]*numColumnPlates))
                plateIDMatrix = np.reshape(self.infoDf['plateName'].values,(plateDimensions[0]*numRowPlates,plateDimensions[1]*numColumnPlates))
                
                blankMatrix = np.reshape(self.allLayouts[0]['blank'].values,(plateDimensions[0]*numRowPlates,plateDimensions[1]*numColumnPlates))

                finalLayoutDict['plateID'] = plateIDMatrix
                finalLayoutDict['wellID'] = wellIDMatrix
                finalLayoutDict['blank'] = blankMatrix
                finalLayoutDict['keys'] = levelValueDict
                
                keyList = []
                for key in finalLayoutDict['keys']:
                    keyList.append(finalLayoutDict['keys'][key])

                fullLayout = np.dstack(keyList)

                uniqueKeys = []
                nonUniquePositions = []
                for row in range(fullLayout.shape[0]):
                    for col in range(fullLayout.shape[1]):
                        wellKeys = fullLayout[row,col,:].tolist()
                        if 'blank' not in wellKeys:
                            if wellKeys not in uniqueKeys:
                                uniqueKeys.append(wellKeys)
                            else:
                                plateID = finalLayoutDict['plateID'][row,col]
                                wellID = finalLayoutDict['wellID'][row,col]
                                nonUniquePositions.append(plateID+'/'+wellID)
                
                if len(nonUniquePositions) == 0:
                    for dt in dataTypeList:
                        with open('misc/layoutDict-'+folderName+'-'+dt+'.pkl','wb') as f:
                            pickle.dump(finalLayoutDict,f)
                        for i in finalLayoutDict['keys']:
                            level = levels[i]
                            createLayoutVisual(baseLayoutDf,finalLayoutDict['keys'][i],i,level,levelValues,plateDimensions,numRowPlates,numColumnPlates,','.join(dataTypeList),infoDf,vlinelist,hlinelist)
                    master.switch_frame(backPage,folderName)
                else:
                    tk.messagebox.showinfo("Duplicate well labels", "These wells (plateID/wellID) have duplicate labels. Please correct them and try again:\n\n"+', '.join(nonUniquePositions))
            else:
                print('Not all levels arranged')
        
        selectionWindow = tk.Frame(self)
        selectionWindow.grid(row=3,column=0)
        advanceVar = tk.StringVar()
        advanceVar.set('yes')
        tk.Button(selectionWindow, text="Deselect Wells",command=lambda: selectWells(False,advanceVar.get())).grid(row=0,column=0)
        tk.Button(selectionWindow, text="Select Wells",command=lambda: selectWells(True,advanceVar.get())).grid(row=0,column=1)
        
        tk.Checkbutton(self,text="Immediately advance level value",variable=advanceVar,onvalue='yes',offvalue='no').grid(row=4,column=0)
        
        editWindow = tk.Frame(self)
        editWindow.grid(row=5,column=0)
        tk.Button(editWindow, text="Copy selection",command=lambda: copySelection()).grid(row=0,column=0)
        tk.Button(editWindow, text="Paste selection",command=lambda: pasteSelection()).grid(row=0,column=1)
        
        actionWindow = tk.Frame(self)
        actionWindow.grid(row=6,column=0)

        actionVar = tk.StringVar()
        actionVar.set('tile')
        tileButton = tk.Radiobutton(actionWindow,text="Tile", variable=actionVar, value='tile')
        repeatButton = tk.Radiobutton(actionWindow,text="Repeat", variable=actionVar, value='repeat')
        tileButton.grid(row=0,column=0,sticky=tk.W)
        repeatButton.grid(row=1,column=0,sticky=tk.W)

        areaVar = tk.StringVar()
        areaVar.set('experiment')
        expButton = tk.Radiobutton(actionWindow,text="Experiment", variable=areaVar, value='experiment')
        plateButton = tk.Radiobutton(actionWindow,text="Plate", variable=areaVar, value='plate')
        expButton.grid(row=0,column=1,sticky=tk.W)
        plateButton.grid(row=1,column=1,sticky=tk.W)
        
        directionVar = tk.StringVar()
        directionVar.set('v')
        vButton = tk.Radiobutton(actionWindow,text="Vertically", variable=directionVar, value='v')
        hButton = tk.Radiobutton(actionWindow,text="Horizontally", variable=directionVar, value='h')
        vButton.grid(row=0,column=2,sticky=tk.W)
        hButton.grid(row=1,column=2,sticky=tk.W)

        tk.Button(actionWindow, text="Perform Action",command=lambda: operateOnSelection(actionVar.get(),areaVar.get(),directionVar.get())).grid(row=0,column=3,rowspan=2)

        tk.Button(self, text="Clear All Wells",command=lambda: clearWells()).grid(row=7,column=0)
        
        buttonWindow = tk.Frame(self)
        buttonWindow.grid(row=8,column=0)
        
        self.FinishButton = tk.Button(buttonWindow, text="Finish",command=lambda: collectInputs(),font='Helvetica 14 bold')
        self.FinishButton.grid(row=0,column=0)
        self.FinishButton['state'] = 'disabled'
        #(self, master,folderName,levels,levelValues,maxNumLevelValues,numRowPlates,numColumnPlates,plateDimensions,dataType,shp)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(BlankSelectionPage,folderName,levels,levelValues,maxNumLevelValues,numColumnPlates,plateDimensions,dataTypeList,secondaryhomepage,backPage)).grid(row=0,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=0,column=2)
