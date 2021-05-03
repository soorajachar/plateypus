#! /usr/bin/env python3
import json,pickle,math,matplotlib,sys,os,string
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import tkinter
import numpy as np
import pandas as pd
sys.path.insert(0, '../plotting/')
import interactiveGUIElements as ipe
from automateCBA import importFCS,beadGate 
import seaborn as sns
from matplotlib import pyplot as plt

cbaGatingParameterDict = {'KDExtrema':['cutoff'],'HDBscan':['min_cluster_size','min_samples']}
cbaGatingParameterBoundsDict = {'cutoff':[0.01,0.1,0.01,0.04],'min_cluster_size':[0,50000,1000,5000],'min_samples':[0,500,10,50]}

class AutomateCBAStartPage(tk.Frame):
    def __init__(self, master,fName,xpNum,x_data,shp,bPage):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        global secondaryhomepage,folderName,expNum,ex_data,backPage
        secondaryhomepage = shp
        folderName = fName
        expNum = xpNum
        ex_data = x_data
        backPage = bPage

        v = tk.StringVar(value='ue')
        rb1a = tk.Radiobutton(mainWindow, text="Add new multiplex gate matrix with dimensions: ",padx = 20, variable=v, value='an')
        t1 = tk.Entry(mainWindow,width=5)
        rb1b = tk.Radiobutton(mainWindow, text="Remove existing multiplex gate matrix",padx = 20, variable=v, value='re')
        rb1c = tk.Radiobutton(mainWindow, text="Use existing multiplex gate matrix",padx = 20, variable=v, value='ue')
        rb1a.grid(row=1,column=0,sticky=tk.W)
        t1.grid(row=1,column=1,sticky=tk.W)
        rb1b.grid(row=2,column=0,sticky=tk.W)
        rb1c.grid(row=3,column=0,sticky=tk.W)
        
        if 'multiplexDict.pkl' not in os.listdir(master.homedirectory+'misc'):
            self.multiplexDict = {}
        else:
            self.multiplexDict = pickle.load(open(master.homedirectory+'misc/multiplexDict.pkl','rb'))
        multiplexes = list(self.multiplexDict.keys())
        self.multiplexMenu = tkinter.ttk.Combobox(mainWindow,values = multiplexes)
        if len(self.multiplexDict) > 0:
            self.multiplexMenu['width'] = len(max(multiplexes,key=len))
        tk.Label(mainWindow,text='Multiplex Sequence: ').grid(row=5,column=0)
        self.multiplexMenu.grid(row=6,column=0)
        
        def collectInput():
            action = v.get()
            if action == 'an':
                matrixDimensions = list(map(int,t1.get().split(',')))
                master.switch_frame(MultiplexCreationPage,matrixDimensions)
            elif action == 're':
                multiplexSequence = self.multiplexMenu.get()
                multiplexDeletion = tk.messagebox.askquestion("Confirmation", "Are you sure you want to delete the multiplex matrix "+multiplexSequence+'?')
                if multiplexDeletion == 'yes':
                    if 'multiplexDict.pkl' in os.listdir(master.homedirectory+'misc'):
                        self.multiplexDict = pickle.load(open(master.homedirectory+'misc/multiplexDict.pkl','rb'))
                        self.multiplexDict = {key:value for key,value in zip(self.multiplexDict.keys(),self.multiplexDict.values()) if key != multiplexSequence}
                        with open(master.homedirectory+'misc/multiplexDict.pkl','wb') as f:
                            pickle.dump(self.multiplexDict,f)
            elif action == 'ue':
                multiplexSequence = self.multiplexMenu.get() 
                multiplexMatrix = pickle.load(open(master.homedirectory+'misc/multiplexDict.pkl','rb'))[multiplexSequence]
                #print(multiplexMatrix)
                tempChannelDict = {'PE ( 561 )-A':'Cytokine Level','APC-A':'Cytokine Identity','BV421-A':'Multiplex 2','FITC-A':'Barcode 1'}
                rawDf = importFCS(tempChannelDict=tempChannelDict)
                print('FCS Files Imported')
                master.switch_frame(MultiplexGatingPage,multiplexMatrix,'KDExtrema',rawDf)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,padx=10,pady=(50,10))
        tk.Button(buttonWindow, text="OK",command=lambda: collectInput()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(secondaryhomepage,folderName,expNum,ex_data,backPage)).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)

class MultiplexCreationPage(tk.Frame):
    def __init__(self, master, matrixDimensions):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10,pady=10)
        
        multiplexFrame = tk.Frame(mainWindow)
        multiplexFrame.pack(side=tk.TOP)
        multiplexEntryList = []
        if len(matrixDimensions) == 1:
            for i in range(matrixDimensions[0]):
                m = tk.Entry(multiplexFrame,width=6)
                m.insert(0, '')
                m.grid(row=0,column=i)
                multiplexEntryList.append(m)
            tk.Label(multiplexFrame,text='Multiplex 1').grid(row=1,column=0,columnspan=matrixDimensions[0])
        elif len(matrixDimensions) == 2:        
            for i in range(matrixDimensions[0]):
                for j in range(matrixDimensions[1]):
                    m = tk.Entry(multiplexFrame,width=6)
                    m.insert(0, '')
                    m.grid(row=i,column=j+1)
                    multiplexEntryList.append(m)
            tk.Label(multiplexFrame,text='Multiplex 1').grid(row=matrixDimensions[0],column=1,columnspan=matrixDimensions[1])
            tk.Label(multiplexFrame,text='Multiplex 2').grid(row=0,column=0,rowspan=matrixDimensions[0])
        else:
            print('This multiplex matrix dimensionality not supported at this time.')
            sys.exit(0)
        
         
        titleWindow = tk.Frame(self)
        titleWindow.pack(side=tk.TOP)
        tk.Label(titleWindow,text='Multiplex Matrix Name: ').grid(row=0,column=0)
        titleEntry = tk.Entry(titleWindow,width=9)
        titleEntry.grid(row=0,column=1)
        
        def collectInput():
            multiplexTitle = titleEntry.get()
            multiplexList = []
            for entry in multiplexEntryList:
                multiplexList.append(entry.get())
            multiplexMatrix = np.array(multiplexList).reshape(matrixDimensions)
            if 'multiplexDict.pkl' not in os.listdir(master.homedirectory+'misc'):
                self.multiplexDict = {multiplexTitle:multiplexMatrix}
            else:
                self.multiplexDict = pickle.load(open(master.homedirectory+'misc/multiplexDict.pkl','rb'))
                self.multiplexDict[multiplexTitle] = multiplexMatrix
            with open(master.homedirectory+'misc/multiplexDict.pkl','wb') as f:
                pickle.dump(self.multiplexDict,f)
            master.switch_frame(AutomateCBAStartPage,folderName,expNum,ex_data,secondaryhomepage,backPage)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP,padx=10,pady=(50,10))
        tk.Button(buttonWindow, text="OK",command=lambda: collectInput()).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(AutomateCBAStartPage,folderName,expNum,ex_data,secondaryhomepage,backPage)).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)

class MultiplexGatingPage(tk.Frame):
    def __init__(self, master, multiplexingMatrix, clusteringMethod, rawData):
        tk.Frame.__init__(self, master)
        
        #Initialize 2x1 canvas for interactive plots
        plotFrame = tk.Frame(self)
        plotFrame.grid(row=0,column=0,columnspan=2)
        fig = plt.figure(figsize=(6, 6))
        gs = fig.add_gridspec(1, 1)
        beadGatingAxis = fig.add_subplot(gs[0])
        self.canvas = FigureCanvasTkAgg(fig,master=plotFrame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()
        
        #Dimensional reduction plot (can be colored/resized/restyled by different level values in dropdowns)
        beadGatingWindow = tk.Frame(self)
        beadGatingWindow.grid(row=1,column=0,sticky=tk.N)
        sliderList = ipe.createParameterAdjustmentSliders(beadGatingWindow,cbaGatingParameterDict[clusteringMethod],cbaGatingParameterBoundsDict)
        
        def updateClusterPlot(sliders):
            beadGatingAxis.clear()
            parametersForClusteringFunction = ipe.getSliderValues(sliderList,cbaGatingParameterDict[clusteringMethod])
            self.clusterdf = beadGate(rawData,**parametersForClusteringFunction)
            clusterPalette = ['blue','grey'] 
            g1 = sns.scatterplot(data=self.clusterdf,x='FSC',y='SSC',s=1,ax=beadGatingAxis,alpha=0.7,hue='Beads',hue_order=['Yes','No'],palette=clusterPalette)
            #beadGatingAxis.set_xlim(self.currentxlims)
            #beadGatingAxis.set_ylim(self.currentylims)
            self.canvas.draw()

        updateClusterPlot(sliderList)
        tk.Button(beadGatingWindow, text="Update bead gate",command=lambda: updateClusterPlot(sliderList)).grid(row=2,column=0)
        
        def exportDataFrames():
            print('Clustered Data Frame And Phenotype Plot Saved')
        
        def okCommand():
            pass
            #exportDataFrames()
            #master.switch_frame(backpage,folderName,secondaryhomepage)

        #Default save and quit buttons
        buttonWindow = tk.Frame(self)
        buttonWindow.grid(row=2,column=0,columnspan=2)
        #tk.Button(buttonWindow, text="OK",command=lambda: master.switch_frame(backpage,folderName,secondaryhomepage)).grid(row=0,column=0)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(ClusteringHomePage,folderName,backpage,secondaryhomepage)).grid(row=0,column=1)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).grid(row=0,column=2)
