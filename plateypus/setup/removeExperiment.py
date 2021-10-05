#!/usr/bin/env python3
import os,subprocess
import tkinter as tk

class RemoveExperimentWindow(tk.Frame):
    def __init__(self,master,backpage):
        tk.Frame.__init__(self, master)
        mainWindow = tk.Frame(self)
        mainWindow.pack(side=tk.TOP,padx=10)
        
        l1 = tk.Label(mainWindow, text="Experiment Name: ")
        experiments = []
        for experimentName in os.listdir('experiments'):
            if '.DS' not in experimentName:
                experiments.append(experimentName)

        selectionVar = tk.StringVar()
        selectionMenu = tk.OptionMenu(mainWindow,selectionVar,*experiments)
        selectionMenu.grid(row=0,column=1)

        def removeExperiment():
            experimentName = selectionMenu.get()
            subprocess.run(['rm','-r','experiments/'+experimentName])

        b = tk.Button(mainWindow,text='Remove experiment',command=lambda:removeExperiment())
        
        l1.grid(row=0,column=0,sticky=tk.W)
        b.grid(row=1,column=0,columnspan=2)

        buttonWindow = tk.Frame(self)
        buttonWindow.pack(side=tk.TOP)
        tk.Button(buttonWindow, text="Back",command=lambda: master.switch_frame(backpage)).pack(side=tk.LEFT)
        tk.Button(buttonWindow, text="Quit",command=lambda: quit()).pack(side=tk.LEFT)
