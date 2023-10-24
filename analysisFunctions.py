import tkinter as tk
import pasdatalib as pdl
import pasphysics as ph
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from numpy import sqrt, array, linspace
from tkinter.filedialog import asksaveasfilename
from scipy.stats import linregress


# Current as of July 5, 2021 at 9:30 PM - Addison Ballif

# This file contains all the analysis and plotting GUI for the PAPPy program. After the parameters
# are collected, we can do additional analysis on them, such as create a S vs W plot. All the analysis
# gui will be in a seperate window, and those seperate windows are contained in this file. 
# 
# The call for the windows to be created is done in the PAPPyIget.py file. 
# 
# 
# The physics for the analysis functions is in the functions. The physics file only contains code for the calculation of the parameters.  

header_font = ('Arial', 18, 'bold') #the font size of the section title
subheader_font = ('Arial', 16)


############# THINGS USEFUL FOR ALL ANALYSIS FUNCTIONS ###########

# ------------- Index Selector --------------
# Was designed for the SW plotter, but could be used for other things as well. 
# 
# It gives you a checklist to select which samples you want to be plotted. 

class IndexSelectorWindow(tk.Frame):
    def __init__(self, parent, master, update_method, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.master = master
        self.parent = parent 

        tk.Label(self, text='Select Samples to Be Plotted: ', font=subheader_font).grid(row=0, column=0)

        self.list_buttons = []
        self.list_variables = []
        self.update_method = update_method

    def clear_list(self):
        #clears all list items
        for item in self.list_buttons:
            item.destroy()

    def populate_list(self):
        #gets data from master and populates list
        self.clear_list()

        for i in range(len(self.master.names)):
            var = tk.IntVar()
            new_button = tk.Checkbutton(self, text=self.master.names[i], command=self.update_method, variable=var)
            new_button.select()
            new_button.grid(row=1+i, column=0)

            self.list_variables.append(var)
            self.list_buttons.append(new_button)

    def get_selected(self):
        #returns a list of indices
        indices = []

        for i in range(len(self.master.names)):
            if self.list_variables[i].get():
                indices.append(i)

        return indices
        

# ----------- Get Matplotlib Savefig filepath function ---------------
# This function's purpose is to get the filepath that we can save a matplotlib figure to. 
# 
def getSaveFilepath():

    files = [('All Files', '*.*'), 
             ('PNG Image', '*.png'),
             ('PDF Document', '*.pdf'),
             ('JPG Image', '*.jpg')]
    return asksaveasfilename(filetypes = files, defaultextension = files)



############################ SW PLOTTER #######################

# ----------- SW Plotter Graph --------
# The graph for the S W Plotter window 
class SWGraph(tk.Frame):
    def __init__(self, parent, master, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        #make the gui here

        #make the graph
        self.figure = plt.Figure(figsize=(7,5), dpi=60)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, self) #add the figure to the window
        self.canvas.get_tk_widget().grid(row=0, column=0)

        self.plotParameters()


    def plotParameters(self):
        #plots s vs w parameters
        self.ax.clear()

        all_x = []
        all_y = []

        if self.parent.control.selectedWParameter.get()=='Right':
            #print('Using Right W Parameter')
            w_list = self.master.rw_param
            w_un_list = self.master.rw_un
        elif self.parent.control.selectedWParameter.get()=='Left':
            #print('Using Left W Parameter')
            w_list = self.master.lw_param
            w_un_list = self.master.lw_un
        elif self.parent.control.selectedWParameter.get()=='Average':
            #assuming it's both (average)
            #print('Using Both W Parameters')
            w_list = (array(self.master.rw_param) + array(self.master.lw_param))/2
            w_un_list = sqrt(array(self.master.rw_un)**2 + array(self.master.lw_un)**2)/2 #independant uncertainty addition
        else:
            #assuming it's both (sum)
            #print('Using Both W Parameters')
            w_list = array(self.master.rw_param) + array(self.master.lw_param)
            w_un_list = sqrt(array(self.master.rw_un)**2 + array(self.master.lw_un)**2) #independant uncertainty addition

        for i in self.parent.control.data_selector.get_selected():
            if self.parent.control.selectedComparisonType.get()=='Parameter Ratio':
                w_un = sqrt((w_un_list[i]/w_list[self.master.ref_idx])**2
                             + (w_list[i]*w_un_list[self.master.ref_idx]/w_list[self.master.ref_idx]**2)**2)
                s_un = sqrt((self.master.s_un[i]/self.master.s_param[self.master.ref_idx])**2
                             + (self.master.s_param[i]*self.master.s_un[self.master.ref_idx]/self.master.s_param[self.master.ref_idx]**2)**2)
                self.ax.errorbar(w_list[i]/w_list[self.master.ref_idx], 
                                 self.master.s_param[i]/self.master.s_param[self.master.ref_idx], 
                                 xerr=w_un, yerr=s_un, 
                                 fmt='.', label=self.master.names[i])

                all_x.append(w_list[i]/w_list[self.master.ref_idx])
                all_y.append(self.master.s_param[i]/self.master.s_param[self.master.ref_idx])
            elif self.parent.control.selectedComparisonType.get()=='Parameter Difference':
                w_un = sqrt(w_un_list[i]**2 + w_un_list[self.master.ref_idx]**2)
                s_un = sqrt(self.master.s_un[i]**2 + self.master.s_un[self.master.ref_idx]**2)
                self.ax.errorbar(array(w_list[i])-w_list[self.master.ref_idx], 
                                 array(self.master.s_param[i])-self.master.s_param[self.master.ref_idx], 
                                 xerr=w_un, yerr=s_un, 
                                 fmt='.', label=self.master.names[i])
                all_x.append(w_list[i]-w_list[self.master.ref_idx])
                all_y.append(self.master.s_param[i]-self.master.s_param[self.master.ref_idx])
            else:
                self.ax.errorbar(w_list[i], self.master.s_param[i], xerr=w_un_list[i], yerr=self.master.s_un[i], 
                                 fmt='.', label=self.master.names[i])

                all_x.append(w_list[i])
                all_y.append(self.master.s_param[i])
            
            

        
        #linear regression
        if self.parent.control.isUsingLinRegress and len(all_x) > 0:
            #linear regression doesn't use error bars. This is because scipy linear regression doesn't have an option for this.
            # Maybe someone in the future could add this.  
            regress = linregress(all_x, all_y)

            x = linspace(min(all_x)-0.1*(max(all_x)-min(all_x)), max(all_x)+0.1*(max(all_x)-min(all_x)), 10)

            self.ax.plot(x, regress.slope*x + regress.intercept, color='black')
            self.ax.text(min(all_x), min(all_y), f' y = ({regress.intercept:.4f} +/- {regress.intercept_stderr:.3f}) + ({regress.slope:.4f} +/- {regress.stderr:.3f}) x')
        elif self.parent.control.isUsingLinRegress and len(all_x) < 1:
            print('Not enough data to run linear regression. ')

        self.ax.set_title('S vs W Plot')
        self.ax.set_xlabel(f'({self.parent.control.selectedWParameter.get()}) W-Parameter')
        self.ax.set_ylabel('S-Parameter')

        if self.parent.control.isShowingLegend:
            self.ax.legend()

        self.canvas.draw()

    def saveImage(self):

        self.figure.savefig(getSaveFilepath(), dpi=300)

# ------------ SW Control
class SWControlWindow(tk.Frame):
    def __init__(self, parent, master, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        tk.Label(self, text='S vs W Plot Settings', font=header_font).grid(row=0, column=0, columnspan=2)

        tk.Label(self, text='Comparison Type: ').grid(row=1, column=0)
        self.selectedComparisonType = tk.StringVar(self)
        self.selectedComparisonType.set('Absolute Parameter')
        self.options_list = ['Absolute Parameter', 'Parameter Difference', 'Parameter Ratio']
        self.optionMenu = tk.OptionMenu(self, self.selectedComparisonType, *self.options_list, command=self.setComparisonType)
        self.optionMenu.grid(row=1, column=1)

        tk.Label(self, text='W Parameter To Use: ').grid(row=2, column=0)
        self.selectedWParameter = tk.StringVar(self)
        self.w_options = ['Right', 'Left', 'Average', 'Sum']
        self.selectedWParameter.set(self.w_options[0])
        self.wOptionsMenu = tk.OptionMenu(self, self.selectedWParameter, *self.w_options, command=self.setSelectedW)
        self.wOptionsMenu.grid(row=2, column=1)

        self.data_selector = IndexSelectorWindow(self, master, update_method=self.updatePlot)
        self.data_selector.grid(row=3, column=0, columnspan=2)
        self.data_selector.populate_list()

        tk.Label(self, text='Run Linear Regression: ').grid(row=4, column=0)
        self.isUsingLinRegress = False
        self.toggleLinRegress = tk.Button(self, text='Turn On', command=self.setIsUsingLinRegress)
        self.toggleLinRegress.grid(row=4, column=1)

        tk.Label(self, text='Plot Legend: ').grid(row=5, column=0)
        self.isShowingLegend = True
        self.toggleShowingLegend = tk.Button(self, text='Turn Off', command=self.setIsShowingLegend)
        self.toggleShowingLegend.grid(row=5, column=1)

        #output file
        tk.Button(self, text='Save Figure', command=self.saveImage).grid(row=6, column=0, columnspan=2)

    def updatePlot(self):
        self.parent.graphWindow.plotParameters()

    def setIsUsingLinRegress(self, **kwargs):
        self.isUsingLinRegress = kwargs.get('value', not self.isUsingLinRegress)
        if self.isUsingLinRegress:
            self.toggleLinRegress.config(text='Turn Off')
        else:
            self.toggleLinRegress.config(text='Turn On')

        self.parent.graphWindow.plotParameters()

    def setIsShowingLegend(self, **kwargs):
        self.isShowingLegend = kwargs.get('value', not self.isShowingLegend)
        if self.isShowingLegend:
            self.toggleShowingLegend.config(text='Turn Off')
        else:
            self.toggleShowingLegend.config(text='Turn On')
        self.parent.graphWindow.plotParameters()

    def setSelectedW(self, selectedParameter):
        #the value is already saved in tkVar
        self.parent.graphWindow.plotParameters()
    
    def setComparisonType(self, selectedType):
        #the type is already updated in the tkVar, so we basically just have to update the plot
        self.parent.graphWindow.plotParameters()

    def saveImage(self):
        #I have to make a new function because when the button is made, the graph has not been created yet. 
        self.parent.graphWindow.saveImage()

        

# ----------- S W Plotter MAIN ---------------
# This is the master window for created SW plots 
class SWPlotterWindow(tk.Toplevel):
    def __init__(self, master, *args, **kwargs):
        tk.Toplevel.__init__(self, *args, **kwargs)
        self.master = master

        ## The SW Plot GUI
        self.title('PAPPy - 2021 - S vs W Plotter')


        self.control = SWControlWindow(self, self.master)
        self.control.grid(row=0, column=1)

        self.graphWindow = SWGraph(self, master)
        self.graphWindow.grid(row=0, column=0)



        self.mainloop() #do this last



############################## DATA PLOTTER ############################### 
# The purpose of the data plotter is to provide a more polished graph of the data. 
# (Mostly I added this to figure out how to add a plotting window, before I added the SW plotter) 

parameter_bounds_alpha = 0.2 #the alpha of the parameter bound fill color

# ----------- GRAPH FRAME ------------
# I'm reusing the graph frame from the main window. It has been copy pasted, and is mostly the same, but has a few adjustments
class DataGraphWindow(tk.Frame):
    def __init__(self, parent, master, dpi=60, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        #make the graph
        self.figure = plt.Figure(figsize=(8,7), dpi=dpi)
        self.figure.subplots_adjust(bottom=0.3, top=0.9)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, self) #add the figure to the window
        self.canvas.get_tk_widget().grid(row=0, column=0)


        #this panel contains a few options
        self.optionsPanel = tk.Frame(self)
        tk.Label(self.optionsPanel, text='Show Background (Hide Other Samples): ').grid(row=0, column=0)
        self.isShowingBackground = True
        self.backgroundToggle = tk.Button(self.optionsPanel, text='Turn Off', command=self.setIsShowingBackground)
        self.backgroundToggle.grid(row=0, column=1)
        

        tk.Label(self.optionsPanel, text='X Units: ').grid(row=1, column=0)
        self.isUsingChannels = False
        self.xUnitToggle = tk.Button(self.optionsPanel, text='Show Channels', command=self.setIsUsingChannels)
        self.xUnitToggle.grid(row=1, column=1)

        tk.Label(self.optionsPanel, text='Show Error Bars: ').grid(row=1, column=2)
        self.isShowingErrorBars = False
        self.errbarToggle = tk.Button(self.optionsPanel, text='Show Error Bars', command=self.setShowErrBars)
        self.errbarToggle.grid(row=1, column=3)

        tk.Label(self.optionsPanel, text='Use Log Scale: ').grid(row=0, column=2)
        self.isUsingLogScale = False
        self.logScaleToggle = tk.Button(self.optionsPanel, text='Turn On', command=self.setUseLogScale)
        self.logScaleToggle.grid(row=0, column=3)

        tk.Label(self.optionsPanel, text='Show Bounds in Legend: ').grid(row=2, column=2)
        self.isShowingBoundsInLegend = True
        self.showBoundsToggle = tk.Button(self.optionsPanel, text='Turn Off', command = self.setIsShowingBoundsInLegend)
        self.showBoundsToggle.grid(row=2, column=3)

        tk.Label(self.optionsPanel, text='Show Integration Regions: ').grid(row=2, column=0)
        self.isShowingRegions = True
        self.showRegionsToggle = tk.Button(self.optionsPanel, text='Turn Off', command = self.setIsShowingRegions)
        self.showRegionsToggle.grid(row=2, column=1)

        tk.Label(self.optionsPanel, text='Zoom in on Peak: ').grid(row=3, column=0)
        self.isZoomedOnCurve = False
        self.zoomCurveToggle = tk.Button(self.optionsPanel, text='Turn On', command=self.setIsZoomedOnCurve)
        self.zoomCurveToggle.grid(row=3, column=1)

        tk.Label(self.optionsPanel, text='Show Grid Lines: ').grid(row=3, column=2)
        self.hasGrid = False
        self.gridToggle = tk.Button(self.optionsPanel, text='Turn On', command=self.setHasGrid)
        self.gridToggle.grid(row=3, column=3)

        self.optionsPanel.grid(row=1, column=0)

        self.plotData()

    def setHasGrid(self, **kwargs):
        self.hasGrid = kwargs.get('value', not self.hasGrid)
        if self.hasGrid:
            self.gridToggle.config(text='Turn Off')
        else:
            self.gridToggle.config(text='Turn On')

        self.plotData()

    def setIsZoomedOnCurve(self, **kwargs):
        self.isZoomedOnCurve = kwargs.get('value', not self.isZoomedOnCurve)
        if self.isZoomedOnCurve:
            self.zoomCurveToggle.config(text='Turn Off')
        else:
            self.zoomCurveToggle.config(text='Turn On')

        self.plotData()

    def setIsShowingBackground(self, **kwargs):
        self.isShowingBackground = kwargs.get('value', not self.isShowingBackground)
        if self.isShowingBackground:
            self.backgroundToggle.config(text='Turn Off')
        else:
            self.backgroundToggle.config(text='Turn On')

        self.plotData()
    def setIsUsingChannels(self, **kwargs):
        #sets the units of the x axis
        self.isUsingChannels = kwargs.get('value', not self.isUsingChannels)
        if self.isUsingChannels:
            self.xUnitToggle.config(text='Show keV')
        else:
            self.xUnitToggle.config(text='Show Channels')

        self.plotData()    
    def setShowErrBars(self, **kwargs):
        #toggles error bars on the graph
        self.isShowingErrorBars = kwargs.get('value', not self.isShowingErrorBars)
        if self.isShowingErrorBars:
            self.errbarToggle.config(text='Hide Error Bars')
        else:
            self.errbarToggle.config(text='Show Error Bars')

        self.plotData()
    def setUseLogScale(self, **kwargs):
        self.isUsingLogScale = kwargs.get('value', not self.isUsingLogScale)
        if self.isUsingLogScale:
            self.logScaleToggle.config(text='Turn Off')
        else:
            self.logScaleToggle.config(text='Turn On')       

        self.plotData()     
    def setIsShowingBoundsInLegend(self, **kwargs):
        self.isShowingBoundsInLegend = kwargs.get('value', not self.isShowingBoundsInLegend)
        if self.isShowingBoundsInLegend:
            self.showBoundsToggle.config(text='Turn Off')
        else:
            self.showBoundsToggle.config(text='Turn On')

        self.plotData()
    def setIsShowingRegions(self, **kwargs):
        self.isShowingRegions = kwargs.get('value', not self.isShowingRegions)
        if self.isShowingRegions:
            self.showRegionsToggle.config(text='Turn Off')
        else:
            self.showRegionsToggle.config(text='Turn On')
        
        self.plotData()


    def plotData(self, **kwargs):

        self.ax.clear()

        if len(kwargs.get('index_list', self.master.data)) < 1:
            #show if no data is loaded
            self.ax.text(0.5, 0.5, 'No Data Loaded', horizontalalignment='center', verticalalignment='center', transform=self.ax.transAxes, bbox=dict(facecolor='red', alpha=0.5))

        #plot the data. updates the graph
        if not self.isUsingChannels:
            
            
            xData = self.master.energy

            if not self.isShowingBackground:
                self.plotAll(xData, index_list=kwargs.get('index_list', range(len(self.master.names))))
            else:
                self.plotReference(xData)

            if self.isShowingRegions:
                #plot a vlines
                lowery = 0
                uppery = self.ax.get_ylim()[1]*2

                self.ax.axvline(self.master.center_keV, lowery, uppery, color='yellow')
                self.ax.vlines([pdl.convertTokeV(self.master.lower, self.master.bins, self.master.energy), 
                                    pdl.convertTokeV(self.master.upper, self.master.bins, self.master.energy)], 
                                    lowery, uppery, color='sienna', label='Curve Bounds' if self.isShowingBoundsInLegend else '_nolegend_')

                self.ax.fill_betweenx([lowery, uppery], pdl.convertTokeV(self.master.parameterWindow.lowerSBound, self.master.bins, self.master.energy), 
                                    pdl.convertTokeV(self.master.parameterWindow.upperSBound, self.master.bins, self.master.energy), 
                                    color='springgreen', alpha=parameter_bounds_alpha, label='S-Param Bounds' if self.isShowingBoundsInLegend else '_nolegend_')
                self.ax.fill_betweenx([lowery, uppery], pdl.convertTokeV(self.master.parameterWindow.lowerLWBound, self.master.bins, self.master.energy), 
                                    pdl.convertTokeV(self.master.parameterWindow.upperLWBound, self.master.bins, self.master.energy), 
                                    color='skyblue', alpha=parameter_bounds_alpha, label='Left W-Param Bounds' if self.isShowingBoundsInLegend else '_nolegend_')
                self.ax.fill_betweenx([lowery, uppery], pdl.convertTokeV(self.master.parameterWindow.lowerRWBound, self.master.bins, self.master.energy), 
                                    pdl.convertTokeV(self.master.parameterWindow.upperRWBound, self.master.bins, self.master.energy), 
                                    color='palevioletred', alpha=parameter_bounds_alpha, label='Right W-Param Bounds' if self.isShowingBoundsInLegend else '_nolegend_')
                
            try:
                if self.isZoomedOnCurve:
                    self.ax.set_xlim(pdl.convertTokeV(self.master.lower, self.master.bins, self.master.energy),
                                    pdl.convertTokeV(self.master.upper, self.master.bins, self.master.energy))
                else:
                    self.ax.set_xlim(min(self.master.energy), max(self.master.energy))
            except ValueError:
                pass #there isn't any data to look at yet. 

            
            self.ax.set_xlabel('Energy (keV)')
            
        else:
            #use channels as the x axis
            
            xData = self.master.bins

            if not self.isShowingBackground:
                self.plotAll(xData, index_list=kwargs.get('index_list', range(len(self.master.names))))
            else:
                self.plotReference(xData)
            
            if self.isShowingRegions:
                #plot a vlines
                lowery = 0
                uppery = self.ax.get_ylim()[1]*2

                self.ax.axvline(self.master.center, lowery, uppery, color='yellow')
                self.ax.vlines([self.master.lower, self.master.upper, ], 
                                    lowery, uppery, color='sienna', label='Curve Bounds' if self.isShowingBoundsInLegend else '_nolegend_')

                self.ax.fill_betweenx([lowery, uppery], self.master.parameterWindow.lowerSBound, self.master.parameterWindow.upperSBound, 
                                    color='springgreen', alpha=parameter_bounds_alpha, label='S-Param Bounds' if self.isShowingBoundsInLegend else '_nolegend_')
                self.ax.fill_betweenx([lowery, uppery], self.master.parameterWindow.lowerLWBound, self.master.parameterWindow.upperLWBound, 
                                    color='skyblue', alpha=parameter_bounds_alpha, label='Left W-Param Bounds' if self.isShowingBoundsInLegend else '_nolegend_')
                self.ax.fill_betweenx([lowery, uppery], self.master.parameterWindow.lowerRWBound, self.master.parameterWindow.upperRWBound, 
                                    color='palevioletred', alpha=parameter_bounds_alpha, label='Right W-Param Bounds' if self.isShowingBoundsInLegend else '_nolegend_')
                

            try:
                if self.isZoomedOnCurve:
                    self.ax.set_xlim(self.master.lower, self.master.upper)
                else:
                    self.ax.set_xlim(min(self.master.bins), max(self.master.bins))
            except ValueError:
                pass #there isn't any data to look at yet. 

            self.ax.set_xlabel('Channels')

            
        #graph settings
        self.ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2)
        self.ax.set_ylabel('Counts')
        if self.isUsingLogScale:
            self.ax.set_yscale('log')
        else:
            self.ax.set_yscale('linear')
        self.canvas.draw()

        self.ax.set_title(self.parent.control.titleEntry.get())


        if self.hasGrid:
            self.ax.grid()
        
        #update the entry fields
        self.parent.control.xlabelEntry.delete(0, tk.END)
        self.parent.control.xlabelEntry.insert(0, self.ax.get_xlabel())

        if len(self.parent.control.ylabelEntry.get()) == 0:
            self.parent.control.ylabelEntry.insert(0, self.ax.get_ylabel())
        else:
            self.ax.set_ylabel(self.parent.control.ylabelEntry.get())

        self.updateYBounds(**kwargs)
        self.updateBoundsLabels()


    def plotAll(self, xData, index_list):
        #use this if you want to plot all files at once. It won't plot backgrounds
        #plot each curve
        for i in index_list:
            if self.isShowingErrorBars:
                self.ax.errorbar(xData, self.master.data[i], label=self.master.names[i], fmt=',', yerr=self.master.data_un[i])
            else:
                self.ax.plot(xData, self.master.data[i], label=self.master.names[i])       
                
    def plotReference(self, xData):
        #use this if you want to plot the reference and background
        
        #plot each curve
        if self.isShowingErrorBars:
            self.ax.errorbar(xData, self.master.data[self.master.ref_idx], label='Reference', yerr=self.master.data_un[self.master.ref_idx], fmt=',')
            self.ax.errorbar(xData, self.master.backgrounds[self.master.ref_idx], label='Background', yerr=self.master.backgrounds_un[self.master.ref_idx])
        else:
            self.ax.plot(xData, self.master.data[self.master.ref_idx], label='Reference')
            self.ax.plot(xData, self.master.backgrounds[self.master.ref_idx], label='Background')

    def updateBoundsLabels(self):
        #updates the Entries in the bounds for the analysis GUI

        self.parent.control.lowerXBound.delete(0, tk.END)
        self.parent.control.upperXBound.delete(0, tk.END)
        self.parent.control.lowerYBound.delete(0, tk.END)
        self.parent.control.upperYBound.delete(0, tk.END)

        self.parent.control.lowerXBound.insert(0, self.ax.get_xlim()[0])
        self.parent.control.upperXBound.insert(0, self.ax.get_xlim()[1])
        self.parent.control.lowerYBound.insert(0, self.ax.get_ylim()[0])
        self.parent.control.upperYBound.insert(0, self.ax.get_ylim()[1])
        

    def updateYBounds(self, **kwargs):
        #updates bounds based off of the reference sample
        max_y = 0
        min_y = 1e19 #an arbitrary high number
        for i in kwargs.get('index_list', range(len(self.master.names))):
            if max(self.master.data[i]) > max_y:
                max_y = max(self.master.data[i])
            if min(self.master.data[i]) < min_y:
                min_y = min(self.master.data[i])

        self.ax.set_ylim(0.8*min_y, 1.2*max_y)
        self.canvas.draw()

    def saveImage(self):

        self.figure.savefig(getSaveFilepath(), dpi=300)

# ------ Control Panel ------------
class DataPlotterControl(tk.Frame):
    def __init__(self, parent, master, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        tk.Label(self, text='Data Plot Settings', font=header_font).grid(row=0, column=0, columnspan=2)

        self.list = IndexSelectorWindow(self, self.master, lambda: self.parent.graph.plotData(index_list=self.list.get_selected()))
        self.list.grid(row=1, column=0, columnspan=2)
        self.list.populate_list()

        tk.Label(self, text='Graph Bounds', font=subheader_font).grid(row=2, column=0, columnspan=2)
        #more controls
        tk.Label(self, text='Lower X Bound: ').grid(row=3, column=0)
        self.lowerXBound = tk.Entry(self)
        self.lowerXBound.grid(row=3, column=1)

        tk.Label(self, text='Upper X Bound: ').grid(row=4, column=0)
        self.upperXBound = tk.Entry(self)
        self.upperXBound.grid(row=4, column=1)

        tk.Label(self, text='Lower Y Bound: ').grid(row=5, column=0)
        self.lowerYBound = tk.Entry(self)
        self.lowerYBound.grid(row=5, column=1)

        tk.Label(self, text='Upper Y Bound: ').grid(row=6, column=0)
        self.upperYBound = tk.Entry(self)
        self.upperYBound.grid(row=6, column=1)

        self.lowerXBound.bind('<Return>', lambda event : self.updateGraphBounds())
        self.upperXBound.bind('<Return>', lambda event : self.updateGraphBounds())
        self.lowerYBound.bind('<Return>', lambda event : self.updateGraphBounds())
        self.upperYBound.bind('<Return>', lambda event : self.updateGraphBounds())

        tk.Label(self, text='Titles and Labels', font=subheader_font).grid(row=7, column=0, columnspan=2)
        
        tk.Label(self, text='Plot Title: ').grid(row=8, column=0)
        self.titleEntry = tk.Entry(self)
        self.titleEntry.grid(row=8, column=1)
        self.titleEntry.bind('<Return>', lambda event : self.setPlotTitle())

        tk.Label(self, text='X Label: ').grid(row=9, column=0)
        self.xlabelEntry = tk.Entry(self)
        self.xlabelEntry.grid(row=9, column=1)
        self.xlabelEntry.bind('<Return>', lambda event : self.setPlotXLabel())
        
        tk.Label(self, text='Y Label: ').grid(row=10, column=0)
        self.ylabelEntry = tk.Entry(self)
        self.ylabelEntry.grid(row=10, column=1)
        self.ylabelEntry.bind('<Return>', lambda event : self.setPlotYLabel())

        #output file
        tk.Button(self, text='Save Figure', command=self.saveImage).grid(row=11, column=0, columnspan=2)


    def setPlotTitle(self):
        self.parent.graph.ax.set_title(self.titleEntry.get())
        self.parent.graph.canvas.draw()
    def setPlotXLabel(self):
        self.parent.graph.ax.set_xlabel(self.xlabelEntry.get())
        self.parent.graph.canvas.draw()
    def setPlotYLabel(self):
        self.parent.graph.ax.set_ylabel(self.ylabelEntry.get())
        self.parent.graph.canvas.draw()

    def updateGraphBounds(self):
        #updates the bounds for the graph
        print('updating bounds')
        #default values if there is no value entered in
        lowerX, upperX = self.parent.graph.ax.get_xlim()
        lowerY, upperY = self.parent.graph.ax.get_ylim()

        try:
            lowerX = float(self.lowerXBound.get())
        except ValueError:
            pass #print("Couldn't Parse Lower X Bound")

        try:
            upperX = float(self.upperXBound.get())
        except ValueError:
            pass #print("Couldn't Parse Upper X Bound")

        try:
            lowerY = float(self.lowerYBound.get())
        except ValueError:
            pass #print("Couldn't Parse Lower Y Bound")

        try: 
            upperY = float(self.upperYBound.get())
        except ValueError:
            pass #print("Couldn't Parse Upper Y Bound")

        self.parent.graph.ax.set_xlim(lowerX, upperX)
        self.parent.graph.ax.set_ylim(lowerY, upperY)

        self.parent.graph.canvas.draw()

    def saveImage(self):
        #I have to make a new function because when the button is made, the graph has not been created yet. 
        self.parent.graph.saveImage()


# ------- THE MAIN WINDOW -----------
class DataPlotterWindow(tk.Toplevel):
    def __init__(self, master, *args, **kwargs):
        tk.Toplevel.__init__(self, *args, **kwargs)
        self.master = master

        self.title('PAPPy 2021 - Data Plotter Tool')

        #make gui here
        self.control = DataPlotterControl(self, self.master)
        self.control.grid(row=0, column=1)

        self.graph = DataGraphWindow(self, self.master, dpi=90) #make this last
        self.graph.grid(row=0, column=0)

        self.mainloop()



# a little something to help people know what to do. 
if __name__ == '__main__':
    print('\nThis is a help file to the PAPPy project. To run the program, run the PAPPy.py file. ')
    