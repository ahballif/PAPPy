import tkinter as tk
import pasdatalib as pdl
import pasphysics as ph
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from os.path import basename
from numpy import array, sqrt
import analysisFunctions as anal #import analysis windows

# Original Code written July 5, 2021 at 9:30 PM - Addison Ballif

# Outline, Plans, and Ideas
# 
# I think for now, I am going to only keep the data that is cropped to the annihilation peak, but
# maybe I should reconsider, because it would be easier to do deconvolution or calibration if we keep
# the other peaks in the data. 
# 
# I am assuming that all the data files have the same calibration, and total number of bins. This way
# I only have to store the x-data once.  

# I am putting window and GUI tools in this file, and physics code in the other file. If it must be related to a GUI,
# then it goes here, but if it is a physics thing, you will find it in one of the other files. The pasdatalib file is
# for importing data and doing data manipulation work, and the pasphysics file contains background subtraction and 
# s-parameter calculation.  

header_font = ('Arial', 18, 'bold') #the font size of the section title
subheader_font = ('Arial', 16)
parameter_bounds_alpha = 0.2 #the alpha of the parameter bound fill color

# ----------- FILE PICKER ------------
# This is the part where you load the csv files. Any number of csv files can be loaded, 
# and the reference sample is also selected.  

class FilePicker(tk.Frame):
    def __init__(self, parent, master, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        self.filepaths = []

        self.statusText = tk.Label(self, text='Please Load a File')
        self.statusText.grid(row=4, column=0, columnspan=2)

        # a dropdown to select the reference file
        tk.Label(self, text='1. Import a Data File', font=header_font).grid(row=0, column=0, columnspan=2)

        #a file browse button
        tk.Label(self, text='Import a CSV File: ').grid(row=1, column=0)
        tk.Button(self, text='Browse', command=self.filebrowse).grid(column=1, row=1)

        #reference sample
        tk.Label(self, text='Select the Reference Sample: ').grid(row=2, column=0)
        self.value_inside = tk.StringVar(self)
        self.value_inside.set('Reference Sample')
        self.optionMenu = tk.OptionMenu(self, self.value_inside, self.master.names, command=self.selectReference)
        self.optionMenu.grid(row=2, column=1)

        #listbox
        self.listbox = tk.Listbox(self, height=3)
        self.listbox.grid(column=0, row=3)
        self.updateListBox()
        tk.Button(self, text='Remove Selected Sample', command=self.removeFile).grid(row=3, column=1)

        
        
        

    def loadFile(self, filepath):
        # This function adds the filepath inputted into the list of filepaths. 
        try: 
            
            #check to see if they already imported it
            if filepath in self.filepaths:
                self.updateStatusText('File Already Loaded. ')
                return

            bins, energy, counts = pdl.importFile(filepath)
            #just override it, because we assume its the same (as in its also starting at 0 and moving up by 1)
            self.master.bins0 = bins 
            #add the energy spectrum to the list of energySpectra
            self.master.energy0.append(energy)
            #add data to the list of data
            self.master.data0.append(counts)
            self.filepaths.append(filepath)
            self.master.names.append(basename(filepath))

            self.updateStatusText('Success! ')
            self.updateOptionsMenu()
            self.updateListBox()
            self.master.updateCroppedData()
            self.master.graph.plotData()
            self.master.bounds.makeGuess()
        except FileNotFoundError:
            self.updateStatusText('File Not Found. ')
            self.filepathEntry.delete(0, 'end')
        except OSError:
            self.updateStatusText('File Not Found. ')
            self.filepathEntry.delete(0, 'end')

    def updateStatusText(self, text):
        #updates the status text
        self.statusText.config(text=text+f'{len(self.master.data0)} files loaded successfully.')
    def updateOptionsMenu(self):
        #makes the option menu
        #updates the selection of options - I am just doing it by deletion and recreation
        self.optionMenu.destroy()
        self.value_inside = tk.StringVar(self)
        try:
            self.value_inside.set(self.master.names[self.master.ref_idx])
        except IndexError:
            self.value_inside.set('Select the Reference File')

        if len(self.master.names) > 0:
            self.optionMenu = tk.OptionMenu(self, self.value_inside, *self.master.names, command=self.selectReference)
        else:
            self.optionMenu = tk.OptionMenu(self, self.value_inside, self.master.names, command=self.selectReference)
        self.optionMenu.grid(row=2, column=1)
    
    def selectReference(self, *args):
        # a function that save the index of the reference sample
        idx = self.master.names.index(self.value_inside.get())
        self.master.ref_idx = idx
        self.master.graph.plotData()

    def filebrowse(self):
        #opens a file brownser
        filepath_list = tk.filedialog.askopenfilenames(initialdir = "",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
        for filepath in filepath_list:
            self.loadFile(filepath=filepath)

    def updateListBox(self):
        #to update the listbox. 
        self.listbox.delete(0,tk.END)#clear 

        for filepath in self.master.names:
            self.listbox.insert(tk.END, filepath)
    def removeFile(self):
        # to remove files which have been loaded. 
        idx = self.listbox.curselection()[0]
        self.master.data0.remove(self.master.data0[idx])
        self.master.energy0.remove(self.master.energy0[idx])
        self.master.names.remove(self.master.names[idx])
        self.filepaths.remove(self.filepaths[idx])
        self.updateListBox()
        self.updateStatusText()
        self.updateOptionsMenu()
        self.master.updateCroppedData()
        self.master.graph.plotData()



# --------------- GRAPH WINDOW ------------------
# The window where we display the data. I think for simplicity right now, 
# I am only going to add the graph which is cropped to  

class GraphWindow(tk.Frame):
    def __init__(self, parent, master, dpi=60, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        #make the graph
        self.figure = plt.Figure(figsize=(8,5.5), dpi=dpi)
        self.figure.subplots_adjust(bottom=0.3, top=0.98)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, self) #add the figure to the window
        self.canvas.get_tk_widget().grid(row=0, column=0)


        #this panel contains a few options
        self.optionsPanel = tk.Frame(self)
        tk.Label(self.optionsPanel, text='Show All Samples: ').grid(row=0, column=0)
        self.isShowingAll = True
        self.allToggle = tk.Button(self.optionsPanel, text='Turn Off', command=self.setIsShowingAll)
        self.allToggle.grid(row=0, column=1)
        

        tk.Label(self.optionsPanel, text='Zoom in on Tails: ').grid(row=0, column=2)
        self.isZoomed = False
        self.zoomToggle = tk.Button(self.optionsPanel, text='Turn On', command=self.setIsZoomed)
        self.zoomToggle.grid(row=0, column=3)

        tk.Label(self.optionsPanel, text='X Units: ').grid(row=1, column=0)
        self.isUsingChannels = False
        self.xUnitToggle = tk.Button(self.optionsPanel, text='Show Channels', command=self.setIsUsingChannels)
        self.xUnitToggle.grid(row=1, column=1)

        tk.Label(self.optionsPanel, text='Show Error Bars: ').grid(row=1, column=2)
        self.isShowingErrorBars = False
        self.errbarToggle = tk.Button(self.optionsPanel, text='Show Error Bars', command=self.setShowErrBars)
        self.errbarToggle.grid(row=1, column=3)

        tk.Label(self.optionsPanel, text='Use Log Scale: ').grid(row=2, column=0)
        self.isUsingLogScale = False
        self.logScaleToggle = tk.Button(self.optionsPanel, text='Turn On', command=self.setUseLogScale)
        self.logScaleToggle.grid(row=2, column=1)

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

    def setIsShowingAll(self, **kwargs):
        self.isShowingAll = kwargs.get('value', not self.isShowingAll)
        if self.isShowingAll:
            self.allToggle.config(text='Turn Off')
        else:
            self.allToggle.config(text='Turn On')

        self.plotData()
    def setIsZoomed(self, **kwargs):
        self.isZoomed = kwargs.get('value', not self.isZoomed)
        if self.isZoomed:
            self.zoomToggle.config(text='Turn Off')
        else:
            self.zoomToggle.config(text='Turn On')
        
        self.updateYBounds()

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

    def plotData(self):

        self.ax.clear()

        if len(self.master.data) < 1:
            #show if no data is loaded
            self.ax.text(0.5, 0.5, 'No Data Loaded', horizontalalignment='center', verticalalignment='center', transform=self.ax.transAxes, bbox=dict(facecolor='red', alpha=0.5))

        #plot the data. updates the graph
        if not self.isUsingChannels and len(self.master.data) > 0:
            
            
            xData = self.master.energy[self.master.ref_idx]

            if self.isShowingAll:
                self.plotAll(xData)
            else:
                self.plotReference(xData)
            
            #plot a vlines
            lowery = 0
            uppery = self.ax.get_ylim()[1]*2

            self.ax.axvline(self.master.center_keV, lowery, uppery, color='yellow')
            self.ax.vlines([pdl.convertTokeV(self.master.lower, self.master.bins, self.master.energy[self.master.ref_idx]), 
                                  pdl.convertTokeV(self.master.upper, self.master.bins, self.master.energy[self.master.ref_idx])], 
                                  lowery, uppery, color='sienna', label='Curve Bounds')

            self.ax.fill_betweenx([lowery, uppery], pdl.convertTokeV(self.master.parameterWindow.lowerSBound, self.master.bins, self.master.energy[self.master.ref_idx]), 
                                  pdl.convertTokeV(self.master.parameterWindow.upperSBound, self.master.bins, self.master.energy[self.master.ref_idx]), 
                                  color='springgreen', alpha=parameter_bounds_alpha, label='S-Param Bounds')
            self.ax.fill_betweenx([lowery, uppery], pdl.convertTokeV(self.master.parameterWindow.lowerLWBound, self.master.bins, self.master.energy[self.master.ref_idx]), 
                                  pdl.convertTokeV(self.master.parameterWindow.upperLWBound, self.master.bins, self.master.energy[self.master.ref_idx]), 
                                  color='skyblue', alpha=parameter_bounds_alpha, label='Left W-Param Bounds')
            self.ax.fill_betweenx([lowery, uppery], pdl.convertTokeV(self.master.parameterWindow.lowerRWBound, self.master.bins, self.master.energy[self.master.ref_idx]), 
                                  pdl.convertTokeV(self.master.parameterWindow.upperRWBound, self.master.bins, self.master.energy[self.master.ref_idx]), 
                                  color='palevioletred', alpha=parameter_bounds_alpha, label='Right W-Param Bounds')
            
            
            #uncomment this if you just want lines rather than the filled region. Gagliardi and Joey Watkins both use a filled region, So I switched to that. 
            # self.ax.vlines([pdl.convertTokeV(self.master.parameterWindow.lowerSBound, self.master.bins, self.master.energy[self.master.ref_idx]), 
            #                       pdl.convertTokeV(self.master.parameterWindow.upperSBound, self.master.bins, self.master.energy[self.master.ref_idx])], 
            #                       lowery, uppery, color='springgreen', label='S-Param Bounds')
            # self.ax.vlines([pdl.convertTokeV(self.master.parameterWindow.lowerLWBound, self.master.bins, self.master.energy[self.master.ref_idx]), 
            #                       pdl.convertTokeV(self.master.parameterWindow.upperLWBound, self.master.bins, self.master.energy[self.master.ref_idx])], 
            #                       lowery, uppery, color='skyblue', label='Left W-Param Bounds')
            # self.ax.vlines([pdl.convertTokeV(self.master.parameterWindow.lowerRWBound, self.master.bins, self.master.energy[self.master.ref_idx]), 
            #                       pdl.convertTokeV(self.master.parameterWindow.upperRWBound, self.master.bins, self.master.energy[self.master.ref_idx])], 
            #                       lowery, uppery, color='palevioletred', label='Right W-Param Bounds')

            #update xbounds to ignore vline positions
            try:
                if self.isZoomedOnCurve:
                    self.ax.set_xlim(pdl.convertTokeV(self.master.lower, self.master.bins, self.master.energy[self.master.ref_idx]),
                                     pdl.convertTokeV(self.master.upper, self.master.bins, self.master.energy[self.master.ref_idx]))
                else:
                    self.ax.set_xlim(min(self.master.energy[self.master.ref_idx]), max(self.master.energy[self.master.ref_idx]))
            except ValueError:
                pass #there isn't any data to look at yet. 

            
            self.ax.set_xlabel('Energy (keV)')
            
        else:
            #use channels as the x axis
            
            xData = self.master.bins

            if self.isShowingAll:
                self.plotAll(xData)
            else:
                self.plotReference(xData)
            
            #plot a vlines
            lowery = 0
            uppery = self.ax.get_ylim()[1]*2

            self.ax.axvline(self.master.center, lowery, uppery, color='yellow')
            self.ax.vlines([self.master.lower, self.master.upper, ], 
                                  lowery, uppery, color='sienna', label='Curve Bounds')

            self.ax.fill_betweenx([lowery, uppery], self.master.parameterWindow.lowerSBound, self.master.parameterWindow.upperSBound, 
                                  color='springgreen', alpha=parameter_bounds_alpha, label='S-Param Bounds')
            self.ax.fill_betweenx([lowery, uppery], self.master.parameterWindow.lowerLWBound, self.master.parameterWindow.upperLWBound, 
                                  color='skyblue', alpha=parameter_bounds_alpha, label='Left W-Param Bounds')
            self.ax.fill_betweenx([lowery, uppery], self.master.parameterWindow.lowerRWBound, self.master.parameterWindow.upperRWBound, 
                                  color='palevioletred', alpha=parameter_bounds_alpha, label='Right W-Param Bounds')
            

            #uncomment this if you just want lines rather than the filled region. Gagliardi and Joey Watkins both use a filled region, So I switched to that. 
            # self.ax.vlines([self.master.parameterWindow.lowerSBound, self.master.parameterWindow.upperSBound], 
            #                       lowery, uppery, color='mediumspringgreen', label='S-Param Bounds')
            # self.ax.vlines([self.master.parameterWindow.lowerLWBound, self.master.parameterWindow.upperLWBound], 
            #                       lowery, uppery, color='skyblue', label='Left W-Param Bounds')
            # self.ax.vlines([self.master.parameterWindow.lowerRWBound, self.master.parameterWindow.upperRWBound], 
            #                       lowery, uppery, color='palevioletred', label='Right W-Param Bounds')

            #update xbounds to ignore vline positions
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

        if self.hasGrid:
            self.ax.grid()


        self.updateYBounds()
        


    def plotAll(self, xData):
        #use this if you want to plot all files at once. It won't plot backgrounds
        #plot each curve
        for i in range(len(self.master.names)):
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

        

    def updateYBounds(self):
        #updates bounds based off of the reference sample
        if not self.isZoomed:
            try:
                self.ax.set_ylim(0, 1.1*max(self.master.data[self.master.ref_idx]))
                self.canvas.draw()
            except IndexError:
                pass #the graph probably just doesn't have anything in it. 
        else:
            try:
                self.ax.set_ylim(0, 0.1*max(self.master.data[self.master.ref_idx]))
                self.canvas.draw()
            except IndexError:
                pass #the graph probably just doesn't have anything in it. 

# -------------- BOUNDS WINDOW ----------------
# this is the window where we set the bounds of the graph, and these are the same bounds 
# use for other calculations
class BoundsWindow(tk.Frame):
    def __init__(self, parent, master, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        #default bounds are alread set in the 

        #bound selection
        tk.Label(self, text='2. Update Curve Bounds', font=header_font).grid(row=0, column=0, columnspan=3)
        
        tk.Label(self, text='Lower Curve Bound (channels): ').grid(row=1, column=0)
        self.lowerEntry = tk.Entry(self)
        self.lowerEntry.grid(row=1, column=1)
        self.lowerEntry.insert(0, self.master.lower)

        tk.Label(self, text='Upper Curve Bound (channels): ').grid(row=2, column=0)
        self.upperEntry = tk.Entry(self)
        self.upperEntry.grid(row=2, column=1)
        self.upperEntry.insert(0, self.master.upper)

        tk.Label(self, text='Background Size (channels): ').grid(row=3, column=0)
        self.background_size_Entry = tk.Entry(self)
        self.background_size_Entry.grid(row=3, column=1)
        self.background_size_Entry.insert(0, self.master.background_size)

        tk.Button(self, text='Update Bounds', command=self.updateXBounds).grid(row=3, column=2)
        tk.Button(self, text='Make a Guess', command=self.makeGuess).grid(row=2, column=2)

        self.lowerEntry.bind('<Return>', lambda event : self.updateXBounds())
        self.upperEntry.bind('<Return>', lambda event : self.updateXBounds())
        self.background_size_Entry.bind('<Return>', lambda event : self.updateXBounds())

        self.centerLabel = tk.Label(self, text=f'Center: {self.master.center_keV:.3f} keV; {self.master.center}th channel')
        self.centerLabel.grid(row=4, column=0)

    def updateGUI (self):
        #updates the labels
        self.centerLabel.config(text=f'Center: {self.master.center_keV:.3f} keV; {self.master.center}th channel')
        self.master.graph.plotData()
        self.master.graph.updateYBounds()

    def updateXBounds(self):
        #updates the bounds
        try:
            lower = int(float(self.lowerEntry.get()))
            upper = int(float(self.upperEntry.get()))
            bgsize = int(float(self.background_size_Entry.get()))

            center = pdl.convertTokeV(float(lower+upper)/2, self.master.bins, self.master.energy[self.master.ref_idx])

            self.master.lower = lower
            self.master.upper = upper
            self.master.center_keV = center
            self.master.center = int(float(upper+lower)/2)
            self.master.background_size = bgsize

            #update the entry so that it can remove decimals
            self.lowerEntry.delete(0, tk.END)
            self.lowerEntry.insert(0, self.master.lower)
            self.upperEntry.delete(0, tk.END)
            self.upperEntry.insert(0, self.master.upper)
            self.background_size_Entry.delete(0, tk.END)
            self.background_size_Entry.insert(0, self.master.background_size)

            self.master.updateCroppedData()
            self.updateGUI()
        except ValueError:
            print('Could Not Parse Bounds. ')

    def makeGuess(self):
        #makes a guess of the bounds
        lower, upper = ph.make_curve_bound_guess(self.master.bins, self.master.energy[self.master.ref_idx], self.master.data[self.master.ref_idx])

        self.master.lower = int(lower)
        self.master.upper = int(upper)
        
        self.lowerEntry.delete(0, tk.END)
        self.upperEntry.delete(0, tk.END)
        self.lowerEntry.insert(0, self.master.lower)
        self.upperEntry.insert(0, self.master.upper)

        center = pdl.convertTokeV(float(lower+upper)/2, self.master.bins, self.master.energy[self.master.ref_idx])
        
        self.master.center_keV = center
        self.master.center = int(float(upper+lower)/2)

        self.master.updateCroppedData()
        self.updateGUI()



# ----------------- BACKGROUND SELECTION WINDOW --------------
class BackgroundSelectionWindow(tk.Frame):
    def __init__(self, parent, master, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.master = parent
        self.master = master

        #make the background gui here
        tk.Label(self, text='3. Background Selection and Subtraction', font=header_font).grid(row=0, column=0, columnspan=2)
        
        tk.Label(self, text='Background Model: ').grid(row=1, column=0)
        self.value_inside = tk.StringVar(self)
        self.value_inside.set(ph.background_models[0])
        self.optionMenu = tk.OptionMenu(self, self.value_inside, *ph.background_models)
        self.optionMenu.grid(row=1, column=1)

        tk.Button(self, text='Calculate Background', command=self.calculateBackground).grid(row=2, column=0)

    def calculateBackground(self):
        # the function that calculates and plots the background.
        for i in range(len(self.master.data)):
            counts = self.master.data[i]
            self.master.backgrounds[i], self.master.backgrounds_un[i] = ph.find_background_shape(self.master.bins, self.master.energy[i], counts, self.master.lower, self.master.upper, self.value_inside.get())
        
        self.master.graph.setIsShowingAll(value=False)
        self.master.background_type = self.value_inside.get()


# --------------- A Table for the S-Paramater Display ---------
class ParamTable(tk.Frame):
    def __init__(self, parent, master, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        self.paramLabelText = 'Parameter'

        #the table will have a name, s and un_s section
        tk.Label(self, text='Sample').grid(row=0, column=0)
        self.namesBox = tk.Listbox(self, height=3)
        self.namesBox.grid(row=1, column=0)

        self.paramLabel = tk.Label(self, text=self.paramLabelText)
        self.paramLabel.grid(row=0, column=1)
        self.paramBox = tk.Listbox(self, height=3)
        self.paramBox.grid(row=1, column=1)
        
        tk.Label(self, text='Uncertainty').grid(row=0, column=2)
        self.unBox = tk.Listbox(self, height=3)
        self.unBox.grid(row=1, column=2)

    def populateTable(self, names, param_list, uncertainty_list):
        #populates the table with the calculated s parameters
        # don't give this function an empty list
        height = len(names)

        #set the heights
        self.namesBox.config(height=height)
        self.paramBox.config(height=height)
        self.unBox.config(height=height)

        #clear the boxes
        self.namesBox.delete(0,tk.END)
        self.paramBox.delete(0,tk.END)
        self.unBox.delete(0,tk.END)

        #add the stuff
        self.namesBox.insert(tk.END, *names)
        self.paramBox.insert(tk.END, *param_list)
        self.unBox.insert(tk.END, *uncertainty_list)

    def set_label(self, paramText):
        self.paramLabelText = paramText
        self.paramLabel.config(text=paramText)

    def clear_table(self):
        #clears the tables and deletes the parameters
        self.namesBox.delete(0,tk.END)
        self.paramBox.delete(0,tk.END)
        self.unBox.delete(0,tk.END)


# ----- Progress Bar for Monte-Carlo Processes --------
# mostly its just a canvas, but I'm adding an update function that updates the canvas 
# according to the self.progress_value 
class ProgressBar(tk.Canvas):
    def __init__(self, parent, *args, **kwargs):
        tk.Canvas.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        self.progress_value = 0 # a float between 0 and 1

        self.n_of_params = 0
        self.completed_params = 0

    def updateProgress(self, **kwargs):
        #udpate the progress
        self.progress_value = self.completed_params/self.n_of_params + kwargs.get('progress_value', 0)/self.n_of_params

        self.create_rectangle(0, 0, self.winfo_width()*self.progress_value, self.winfo_height(), fill='blue')
        self.update()

    def setNumParameters(self, n_of_params):
        self.n_of_params = n_of_params
    def setCompletedParams(self, completed_params):
        self.completed_params = completed_params

    def reset_progress(self):
        #resets the progress graphic
        self.create_rectangle(0, 0, self.winfo_width(), self.winfo_height(), fill='white')
        self.progress_value = 0
        



# -------------------- S-Parameter Selection Window ------------
# The s-parameter selection window contains the bounds for the s, left w, and right w parameters. 
# 
# The GUI is made up of entries for each of these parameters. There is also a button that makes a guess for each parameter. 
# Also contained in this ParameterWindow are settings for the parameter calculation, and the parameter output table. 
# 
# The results from the calculation are placed in the ParameterTable as well as stored in the parameter lists which
# are contained in the MainApplication class. 
# 

 
class ParameterWindow(tk.Frame):
    def __init__(self, parent, master, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        self.lowerSBound = 0
        self.upperSBound = 0

        self.lowerLWBound = 0
        self.upperLWBound = 0
        
        self.lowerRWBound = 0
        self.upperRWBound = 0

        #make the s-parameter gui here. 
        tk.Label(self, text='4. Finding Parameters', font=header_font).grid(row=0, column=0, columnspan=3)
        
        #s-parameter bounds
        #tk.Label(self, text='S-Parameter', font=subheader_font).grid(row=1, column=0, columnspan=3)
        tk.Label(self, text='Lower S Bound (channels): ').grid(row=2, column=0)
        self.lowerSEntry = tk.Entry(self)
        self.lowerSEntry.grid(row=2, column=1)
        tk.Label(self, text='Upper S Bound (channels): ').grid(row=3, column=0)
        self.upperSEntry = tk.Entry(self)
        self.upperSEntry.grid(row=3, column=1)
        tk.Button(self, text='Make a Guess', command=self.makeGuessS).grid(row=3, column=2)
        self.lowerSEntry.bind('<Return>', lambda event: self.updateBounds())
        self.upperSEntry.bind('<Return>', lambda event: self.updateBounds())

        #left w-parameter bounds
        #tk.Label(self, text='Left W-Parameter', font=subheader_font).grid(row=4, column=0, columnspan=3)
        tk.Label(self, text='Lower Left W Bound (channels): ').grid(row=5, column=0)
        self.lowerLWEntry = tk.Entry(self)
        self.lowerLWEntry.grid(row=5, column=1)
        tk.Label(self, text='Upper Left W Bound (channels): ').grid(row=6, column=0)
        self.upperLWEntry = tk.Entry(self)
        self.upperLWEntry.grid(row=6, column=1)
        tk.Button(self, text='Make a Guess', command=self.makeGuessLW).grid(row=6, column=2)
        self.lowerLWEntry.bind('<Return>', lambda event: self.updateBounds())
        self.upperLWEntry.bind('<Return>', lambda event: self.updateBounds())

        #right w-parameter bounds
        #tk.Label(self, text='Right W-Parameter', font=subheader_font).grid(row=7, column=0, columnspan=3)
        tk.Label(self, text='Lower Right W Bound (channels): ').grid(row=8, column=0)
        self.lowerRWEntry = tk.Entry(self)
        self.lowerRWEntry.grid(row=8, column=1)
        tk.Label(self, text='Upper Right W Bound (channels): ').grid(row=9, column=0)
        self.upperRWEntry = tk.Entry(self)
        self.upperRWEntry.grid(row=9, column=1)
        tk.Button(self, text='Make a Guess', command=self.makeGuessRW).grid(row=9, column=2)
        self.lowerRWEntry.bind('<Return>', lambda event: self.updateBounds())
        self.upperRWEntry.bind('<Return>', lambda event: self.updateBounds())

        #some settings
        #tk.Label(self, text='Calculation Settings', font=subheader_font).grid(row=10, column=0, columnspan=3)
        tk.Label(self, text='Use Trapezoidal Integration: ').grid(row=11, column=0)
        self.isUsingTrapz = False
        self.trapzToggle = tk.Button(self, text='Turn On', command=self.setIsUsingTrapz)
        self.trapzToggle.grid(row=11, column=1)

        tk.Label(self, text='Number of Monte-Carlo Iterations: ').grid(row=12, column=0)
        self.iterationsEntry = tk.Entry(self)
        self.iterationsEntry.insert(0, 500)
        self.iterationsEntry.grid(row=12, column=1)

        outPutRow = 13 #the row where output and calculate button starts

        # the BUTTON
        tk.Button(self, text='Calculate Parameters', command=self.calculateSParameters, font=('Corbel', 15, 'bold')).grid(row=outPutRow, column=0, columnspan=3)
        
        self.progressText = tk.Label(self, text='Parameters Calculated: 0/0')
        self.progressText.grid(row=outPutRow+1, column=1)
        self.progressBar = ProgressBar(self, bg='white', height=20, width=200)
        self.progressBar.grid(row=outPutRow+1, column=2)

        #output tables
        #tk.Label(self, text='S-Parameters').grid(row=outPutRow+2, column=0, columnspan=3)
        self.s_table = ParamTable(self, self.master)
        self.s_table.grid(row=outPutRow+3, column=0, columnspan=3)
        self.s_table.set_label('S-Parameter')

        #tk.Label(self, text='Left W-Parameters').grid(row=outPutRow+4, column=0, columnspan=3)
        self.lw_table = ParamTable(self, self.master)
        self.lw_table.grid(row=outPutRow+5, column=0, columnspan=3)
        self.lw_table.set_label('Left W-Parameter')

        #tk.Label(self, text='Right W-Parameters').grid(row=outPutRow+6, column=0, columnspan=3)
        self.rw_table = ParamTable(self, self.master)
        self.rw_table.grid(row=outPutRow+7, column=0, columnspan=3)
        self.rw_table.set_label('Right W-Parameter')
    
    def updateBounds(self):
        #used to update graph and bound values as you change the entry

        #s-parameter
        try:
            self.lowerSBound = int(self.lowerSEntry.get())
            self.upperSBound = int(self.upperSEntry.get())

            self.lowerSEntry.delete(0, 'end')
            self.upperSEntry.delete(0, 'end')
            
            self.lowerSEntry.insert(0, self.lowerSBound)
            self.upperSEntry.insert(0, self.upperSBound)
        except ValueError:
            #print("Couldn't Parse S-Parameter Bound Input")

            self.lowerSBound = 0
            self.upperSBound = 0
            
            self.lowerSEntry.delete(0, 'end')
            self.upperSEntry.delete(0, 'end')

        #left w parameter
        try:
            self.lowerLWBound = int(self.lowerLWEntry.get())
            self.upperLWBound = int(self.upperLWEntry.get())

            self.lowerLWEntry.delete(0, 'end')
            self.upperLWEntry.delete(0, 'end')
            
            self.lowerLWEntry.insert(0, self.lowerLWBound)
            self.upperLWEntry.insert(0, self.upperLWBound)
        except ValueError:
            #print("Couldn't Parse Left W-Parameter Bound Input")

            self.lowerLWBound = 0
            self.upperLWBound = 0

            self.lowerLWEntry.delete(0, 'end')
            self.upperLWEntry.delete(0, 'end')

        #right w parameter
        try:
            self.lowerRWBound = int(self.lowerRWEntry.get())
            self.upperRWBound = int(self.upperRWEntry.get())

            self.lowerRWEntry.delete(0, 'end')
            self.upperRWEntry.delete(0, 'end')
            
            self.lowerRWEntry.insert(0, self.lowerRWBound)
            self.upperRWEntry.insert(0, self.upperRWBound)
        except ValueError:
            #print("Couldn't Parse Right W-Parameter Bound Input")

            self.lowerRWBound = 0
            self.upperRWBound = 0

            self.lowerRWEntry.delete(0, 'end')
            self.upperRWEntry.delete(0, 'end')

        self.master.graph.plotData() #update display

    def setIsUsingTrapz(self, **kwargs):
        self.isUsingTrapz = kwargs.get('value', not self.isUsingTrapz)
        if self.isUsingTrapz:
            self.trapzToggle.config(text='Turn Off')
        else:
            self.trapzToggle.config(text='Turn On')

    def makeGuessS(self):
        #makes a guess for p-value bounds based off of reference sample
        lower, upper = ph.make_s_parameter_bound_guess(self.master.bins, self.master.data[self.master.ref_idx])

        self.lowerSBound = int(lower)
        self.upperSBound = int(upper)

        self.lowerSEntry.delete(0, 'end')
        self.upperSEntry.delete(0, 'end')
        
        self.lowerSEntry.insert(0, int(lower))
        self.upperSEntry.insert(0, int(upper))
        
        self.master.graph.plotData()

    def makeGuessLW(self):
        #makes a guess for p-value bounds based off of reference sample
        lower, upper = ph.make_left_w_parameter_bound_guess(self.master.bins, self.master.data[self.master.ref_idx], (self.master.lower, self.master.upper))

        self.lowerLWBound = int(lower)
        self.upperLWBound = int(upper)

        self.lowerLWEntry.delete(0, 'end')
        self.upperLWEntry.delete(0, 'end')
        
        self.lowerLWEntry.insert(0, int(lower))
        self.upperLWEntry.insert(0, int(upper))
        
        self.master.graph.plotData()
    
    def makeGuessRW(self):
        #makes a guess for p-value bounds based off of reference sample
        lower, upper = ph.make_right_w_parameter_bound_guess(self.master.bins, self.master.data[self.master.ref_idx], (self.master.lower, self.master.upper))

        self.lowerRWBound = int(lower)
        self.upperRWBound = int(upper)

        self.lowerRWEntry.delete(0, 'end')
        self.upperRWEntry.delete(0, 'end')
        
        self.lowerRWEntry.insert(0, int(lower))
        self.upperRWEntry.insert(0, int(upper))
        
        self.master.graph.plotData()

    ### CALCULATING PARAMETERS ###
    def calculateSParameters(self):

        self.master.backgroundWindow.calculateBackground() #calculate the background
        
        #calculates all the s and w parameters
        s_param = []
        s_un = []
        lw_param = []
        lw_un = []
        rw_param = []
        rw_un = []
        names = self.master.names

        self.updateBounds() #this checks for empty fields

        #for status text
        n_of_params = 3
        if self.lowerSBound == 0 and self.upperSBound == 0:
            n_of_params -= 1
        if self.lowerLWBound == 0 and self.upperLWBound == 0:
            n_of_params -= 1
        if self.lowerRWBound == 0 and self.upperRWBound == 0:
            n_of_params -= 1
        n_of_params *= len(names)

        params_done = 0

        self.progressText.config(text=f'Parameters Calculated: {0}/{n_of_params}') #update progress bar
        self.progressBar.setCompletedParams(0)
        self.progressBar.setNumParameters(n_of_params)
        self.progressBar.reset_progress()

        #get all settings and data before looping. That way if it gets changed, it won't mess up the calculation
        curve_bounds = (self.master.lower, self.master.upper)
        s_bounds = (self.lowerSBound, self.upperSBound)
        lw_bounds = (self.lowerLWBound, self.upperLWBound)
        rw_bounds = (self.lowerRWBound, self.upperRWBound)

        background_type = self.master.background_type
        use_trapz = self.isUsingTrapz

        d_bins = self.master.bins
        d_data = self.master.data
        d_background = self.master.backgrounds
        d_background_un = self.master.backgrounds_un

        try:
            n_of_iterations = int(self.iterationsEntry.get())
        except ValueError:
            if self.master.background_type == 'erf':
                print('Could not Parse Number of Iterations, assuming 500. ')
                n_of_iterations = 500

        #loop through data sets
        for i in range(len(names)):
               
            if s_bounds[0] != 0 and s_bounds[1] != 0:
                params_done += 1
                s_parameter, s_parameter_un = ph.calculate_parameter(d_bins, d_data[i], d_background[i], d_background_un[i], 
                                                        curve_bounds, s_bounds, background_type, trapz=use_trapz, progress_bar=self.progressBar, number_of_iterations=n_of_iterations)
                s_param.append(s_parameter)
                s_un.append(s_parameter_un)
            else:
                s_param.append(0)
                s_un.append(0)

            self.progressText.config(text=f'Parameters Calculated: {params_done}/{n_of_params}') #update progress bar
            self.progressBar.setCompletedParams(params_done)

            if lw_bounds[0] != 0 and lw_bounds[1] != 0:
                params_done+=1
                lw_parameter, lw_parameter_un = ph.calculate_parameter(d_bins, d_data[i], d_background[i], d_background_un[i], 
                                                        curve_bounds, lw_bounds, background_type, trapz=use_trapz, progress_bar=self.progressBar, number_of_iterations=n_of_iterations)
                lw_param.append(lw_parameter)
                lw_un.append(lw_parameter_un)
            else:
                lw_param.append(0)
                lw_un.append(0)

            self.progressText.config(text=f'Parameters Calculated: {params_done}/{n_of_params}') #update progress bar
            self.progressBar.setCompletedParams(params_done)

            if rw_bounds[0] != 0 and rw_bounds[1] != 0:
                params_done+=1
                rw_parameter, rw_parameter_un = ph.calculate_parameter(d_bins, d_data[i], d_background[i], d_background_un[i], 
                                                        curve_bounds, rw_bounds, background_type, trapz=use_trapz, progress_bar=self.progressBar, number_of_iterations=n_of_iterations)
                rw_param.append(rw_parameter)
                rw_un.append(rw_parameter_un)
            else:
                rw_param.append(0)
                rw_un.append(0)

            self.progressText.config(text=f'Parameters Calculated: {params_done}/{n_of_params}') #update progress bar
            self.progressBar.setCompletedParams(params_done)



        self.s_table.populateTable(names, s_param, s_un)
        self.lw_table.populateTable(names, lw_param, lw_un)
        self.rw_table.populateTable(names, rw_param, rw_un)

        self.master.s_param = s_param
        self.master.s_un = s_un
        self.master.lw_param = lw_param
        self.master.lw_un = lw_un
        self.master.rw_param = rw_param
        self.master.rw_un = rw_un

# ------------------------ PLOTS AND ADDITIONAL ANALYSIS -------
# Sometimes it is useful to do additional analysis to the parameters inside the program. 
# This window contains buttons which will run various forms of additional analysis. 
# 
# When a button is clicked, it will open a new TopLevel window. All additional analysis GUI will 
# be contained in a new window.  

class AnalysisWindow(tk.Frame):
    def __init__(self, parent, master, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.master = master

        #additional function buttons
        tk.Label(self, text='5. Plots and Analysis', font=header_font).grid(row=0, column=0, columnspan=3)
        
        tk.Button(self, text='S vs W Plot', command=self.s_rw_plot).grid(row=1, column=0)
        tk.Button(self, text='Data Plotter', command=self.data_plot).grid(row=1, column=1)
        tk.Button(self, text='Export Parameters as CSV', command=self.exportParameters).grid(row=1, column=2)
    
        

    def s_rw_plot(self):
        #make an s vs w plot using the right w parameter
        anal.SWPlotterWindow(self.master)

    def data_plot(self):
        #makes a more presentable version of the data plot from the home screen
        anal.DataPlotterWindow(self.master)

    def exportParameters(self):
        #exports the parameters in csv form. 
        files = [('All Files', '*.*'), 
             ('CSV File', '*.csv')]
        filepath = tk.filedialog.asksaveasfilename(filetypes = files, defaultextension = files)
        pdl.exportParameters(filepath, ['Sample Name', 'S Parameter', 'S Uncertainty', 'Left W Parameter', 'Left W Uncertainty', 'Right W Parameter', 'Right W Uncertainty'], 
                                            self.master.names,
                                            self.master.s_param, self.master.s_un, 
                                            self.master.lw_param, self.master.lw_un,
                                            self.master.rw_param, self.master.rw_un)
    
        

         
    
# ------------------------- MAIN APPLICATION --------------------------
# The main window houses the different panels of the GUI. It also is the class
# that stores all the information. 
# 
# The original data is stored in the bins0, energy0, data0 lists. bins and energy are 1d lists with the x 
# data. The data0 list contains each y-data set in list form. It is an NxM matrix, where N is the number of files
# and M is the number of channels in each file. 
# 
# The bins, energy, and data lists are cropped version of the original data. Changing the bounds in the BoundsWindow
# class will re-calculate these lists. There is also background, data_un, background_un and lists. The uncertainty for the
# data is only calculated for the cropped data for simplicity. The background and background_un lists contain the background 
# signal. After choosing the background model and calculating background, it is stored in this NxM matrix. 
# 
# The parameters and uncertainties are stored in a list of size N. This is done after the parameters are calculated, and is mainly
# important for use in the Plots and Analysis section. The updateCropped data, which recalculates the cropped bins, energy, and data 
# lists also resets these parameters lists, as well as the background and uncertainty lists.  
# 
# The background type is important to remember, because it is used in the parameter calculation. This is because different background
# models will require different methods of uncertainty propegation. After the background model is selected and background is calculated,
# the background type will be stored in this variable inside the main application. 
# 
# The curve bounds are stored in this class as well. The center and center_keV are only dispayed and never used. These are to help
# the user identify good bounds. 
# 
# The lower, upper, and background_size variables are edited by the bounds window. Parameter bounds are stored in the parameter window,
# but curve bounds and background size are stored here. 
# 
# The rest of the variables in this class are the different GUI windows, and the tk.Frames that organize the GUI. 
#  
class MainApplication(tk.Frame):
    def __init__(self, parent, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        # These lists store the data from the files
        # There is a 0 version which contains the entire spectrum
        # The bins, energy, data lists constain the cropped data. 
        self.bins0 = [] # a single array
        self.energy0 = [] # this will be an array of arrays of energy spectra (because different samples may have different calibrations)
        self.data0 = [] #this will be an array of arrays of counts
        
        self.bins = []
        self.energy = []
        self.data = []
        self.backgrounds = [] #the background is also an array of arrays
        self.data_un = []
        self.backgrounds_un = []

        self.s_param = []
        self.s_un = []
        self.lw_param = []
        self.lw_un = []
        self.rw_param = []
        self.rw_un = []

        self.ref_idx = 0 #the index of the reference
        self.names = []

        self.background_type = 'none' #the type of background. Used for parameter calculations

        #bounds variables
        self.center_keV = 511 #keV - for reading only, we can't set this in the GUI
        self.center = 5110 #channels - for reading only, we can't set this in the GUI

        self.lower = 5050 #channels
        self.upper = 5165 #channels
        self.background_size = 100 #channels

        #make the gui here. 
        # the left side and right side panels are used so that the rows in each column do not need to line up. 
        self.leftSide = tk.Frame(self)
        self.rightSide = tk.Frame(self)
        self.leftSide.grid(row=0, column=0)
        self.rightSide.grid(row=0, column=1)

        self.filePicker = FilePicker(self.leftSide, self)
        self.filePicker.grid(row=0, column=0)

        self.parameterWindow = ParameterWindow(self.rightSide, self) 
        self.parameterWindow.grid(row=1, column=0)
        
        self.bounds = BoundsWindow(self.leftSide, self)
        self.bounds.grid(row=2, column=0)

        self.backgroundWindow = BackgroundSelectionWindow(self.rightSide, self)
        self.backgroundWindow.grid(row=0, column=0)

        self.analysisWindow = AnalysisWindow(self.rightSide, self)
        self.analysisWindow.grid(row=2, column=0)

            #we must initialize graph last to make sure we have all the information needed
        self.graph = GraphWindow(self.leftSide, self) 
        self.graph.grid(row=1, column=0)


    
    def updateCroppedData(self):
        #updates the cropped data
        self.data = []
        self.energy = []
        self.backgrounds = []
        self.data_un = []
        self.backgrounds_un = []
        self.background_type = 'none'

        self.s_param = []
        self.s_un = []
        self.lw_param = []
        self.lw_un = []
        self.rw_param = []
        self.rw_un = []

        for counts0 in self.data0:
            self.bins, energySpectrum, counts = pdl.cropWindow(self.bins0, self.energy0[self.data0.index(counts0)], counts0, self.lower-self.background_size, self.upper+self.background_size)
            self.data.append(counts)
            self.energy.append(energySpectrum)
            self.data_un.append(sqrt(array(counts)))
            self.backgrounds.append(array(counts)*0)
            self.backgrounds_un.append(array(counts)*0)

            self.s_param.append(0)
            self.lw_param.append(0)
            self.rw_param.append(0)
            self.s_un.append(0)
            self.lw_un.append(0)
            self.rw_un.append(0)

def onFrameConfigure(canvas):
    '''Reset the scroll region to encompass the inner frame'''
    canvas.configure(scrollregion=canvas.bbox("all"))

if __name__ == '__main__':
    root = tk.Tk()

    MAX_WINDOW_WIDTH = 1200
    MAX_WINDOW_HEIGHT = 755

    width  = root.winfo_screenwidth()
    height = root.winfo_screenheight()
    root.geometry(f'{(MAX_WINDOW_WIDTH if width > MAX_WINDOW_WIDTH else width)}x{(MAX_WINDOW_HEIGHT if height > MAX_WINDOW_HEIGHT else height)}')

    #putting the main application in a scrollwindow for smaller screens
    canvas = tk.Canvas(root, borderwidth=0)
    main = MainApplication(canvas)

    vsb = tk.Scrollbar(root, orient='vertical', command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)

    hsb = tk.Scrollbar(root, orient='horizontal', command=canvas.xview)
    canvas.configure(xscrollcommand=hsb.set)

    vsb.pack(side='right', fill='y')
    hsb.pack(side='bottom', fill='x')


    canvas.pack(side='left', fill='both', expand=True)
    canvas.create_window((0,0), window=main, anchor='nw')

    main.bind('<Configure>', lambda event, canvas=canvas : onFrameConfigure(canvas))

    #load some sample data
    main.filePicker.loadFile('example_reference.csv')
    main.filePicker.loadFile('example_test.csv')

    root.title('PAPPy - 2021')
    root.mainloop()

    
    
