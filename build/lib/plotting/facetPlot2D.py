#!/usr/bin/env python3
import seaborn as sns
    
def plot(plottingDf,subsettedDf,kwargs,facetKwargs,auxillaryKwargs,plotOptions):
    #Make sure there are markers at each column variable
    if 'style' not in kwargs.keys():
        fg = sns.relplot(data=plottingDf,marker='o',kind=auxillaryKwargs['subPlotType'],facet_kws=facetKwargs,ci=False,**kwargs,**plotOptions['X']['figureDimensions'],**auxillaryKwargs['cmap'])
    else:
        fg = sns.relplot(data=plottingDf,markers=True,kind=auxillaryKwargs['subPlotType'],facet_kws=facetKwargs,ci=False,**kwargs,**plotOptions['X']['figureDimensions'],**auxillaryKwargs['cmap'])
        
    #X and Y Axis Scaling for 2D plots
    for axis in plotOptions:
        k = len(fg.fig.get_axes())
        if 'Y' in axis:
            if plotOptions[axis]['axisScaling'] == 'Logarithmic':
                for i in range(k):
                    fg.fig.get_axes()[i].set_yscale('log')
            elif plotOptions[axis]['axisScaling'] == 'Biexponential':
                for i in range(k):
                    fg.fig.get_axes()[i].set_yscale('symlog',linthreshx=plotOptions[axis]['linThreshold'])
            
            if str(plotOptions[axis]['limit'][0]) != '' or str(plotOptions[axis]['limit'][1]) != '':
                for i in range(k):
                    if str(plotOptions[axis]['limit'][0]) != '' and str(plotOptions[axis]['limit'][1]) != '':
                        fg.fig.get_axes()[i].set_ylim(bottom=float(plotOptions[axis]['limit'][0]),top=float(plotOptions[axis]['limit'][1]))
                    else:
                        if str(plotOptions[axis]['limit'][0]) != '':
                            fg.fig.get_axes()[i].set_ylim(bottom=float(plotOptions[axis]['limit'][0]))
                        else:
                            fg.fig.get_axes()[i].set_ylim(top=float(plotOptions[axis]['limit'][1]))
        else:
            if plotOptions[axis]['axisScaling'] == 'Logarithmic':
                for i in range(k):
                    fg.fig.get_axes()[i].set_xscale('log')
            elif plotOptions[axis]['axisScaling'] == 'Biexponential':
                for i in range(k):
                    fg.fig.get_axes()[i].set_xscale('symlog',linthreshx=plotOptions[axis]['linThreshold']) 
            
            if str(plotOptions[axis]['limit'][0]) != '' or str(plotOptions[axis]['limit'][1]) != '':
                for i in range(k):
                    if str(plotOptions[axis]['limit'][0]) != '' and str(plotOptions[axis]['limit'][1]) != '':
                        fg.fig.get_axes()[i].set_xlim(bottom=float(plotOptions[axis]['limit'][0]),top=float(plotOptions[axis]['limit'][1]))
                    else:
                        if str(plotOptions[axis]['limit'][0]) != '':
                            fg.fig.get_axes()[i].set_xlim(bottom=float(plotOptions[axis]['limit'][0]))
                        else:
                            fg.fig.get_axes()[i].set_xlim(top=float(plotOptions[axis]['limit'][1])) 
    return fg
