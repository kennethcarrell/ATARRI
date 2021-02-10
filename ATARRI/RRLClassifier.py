import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from lightkurve import search_tesscut
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
import tkinter
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import warnings
warnings.filterwarnings("ignore",category=UserWarning)

class RRLClassifier:
    '''A GUI for visualizing and saving information about an RR Lyrae star from TESS FFI data.'''

    # Initialize variables and GUI
    def __init__(self, master):
        self.RA = []
        self.DE = []
        self.tpf = 0
        self.apthresh = 3.0
        self.target_mask0 = 0
        self.bkgr_mask0 = 0
        self.star_lc = 0
        self.lspg = 0
        self.pwlspg = 0
        self.PERIOD = []
        self.TESSID = []
        self.NSECS = []
        self.RTYPE = []
        self.APTEST = []
        self.APSIZE = []
        self.PHTEST = []
        self.PDR = []
        self.NRR = []
        self.BLA = []
        self.ORA = []
        self.ODE = []
        self.NOTES = []
        self.NAME = []
        self.NDONE = 0
        self.tableOut = Table()
        self.starName = ''
        self.infile = ''
        self.outfile = ''

        # Create Frames
        self.master = master
        self.master.title("ATARRI GUI")
        self.infoframe = tkinter.Frame(self.master)
        self.plotframe = tkinter.Frame(self.master)
        self.optionsframe = tkinter.Frame(self.master)
        self.apselectframe = tkinter.Frame(self.master)
        self.notesframe = tkinter.Frame(self.master)

        # Information areas
        self.fNameButton = tkinter.Button(self.infoframe, text="Open File", fg="blue", command=self.retrieveInput)
        self.fNameButton.grid(row=0,column=0,padx=80)
        self.fNameLabel = tkinter.Label(self.infoframe, text="Input File Name: ")
        self.fNameLabel.grid(row=0,column=1,padx=80)
        self.fNumberLabel = tkinter.Label(self.infoframe, text="N Stars = 0")
        self.fNumberLabel.grid(row=0,column=2,padx=80)
        self.saveButton = tkinter.Button(self.infoframe, text="Save This Star", fg="green", command=self.saveStar)
        self.saveButton.grid(row=0,column=3,padx=80)
        self.inspectLCButton = tkinter.Button(self.infoframe, text="Inspect Lightcurve", fg="blue", command=self.inspectLC)
        self.inspectLCButton.grid(row=0,column=4,padx=80)
        self.sTICLabel = tkinter.Label(self.infoframe, text="TIC: ")
        self.sTICLabel.grid(row=1,column=0)
        self.sPosLabel = tkinter.Label(self.infoframe, text="(RA,DEC): ")
        self.sPosLabel.grid(row=1,column=1)
        self.inspectPhasedButton = tkinter.Button(self.infoframe, text="Inspect Phased Plot", fg="blue", command=self.inspectPhased)
        self.inspectPhasedButton.grid(row=1,column=4)

        # Classification Type
        self.RRLType = tkinter.IntVar()
        self.RRLType.set(0)
        self.typeLabel = tkinter.Label(self.optionsframe, text="Select Type")
        self.typeLabel.grid(row=0,column=0,sticky=tkinter.W,padx=30)
        self.typeabButton = tkinter.Radiobutton(self.optionsframe, text="RRab", variable=self.RRLType, value=1)
        self.typeabButton.grid(row=1,column=0,sticky=tkinter.W,padx=35)
        self.typecButton = tkinter.Radiobutton(self.optionsframe, text="RRc", variable=self.RRLType, value=2)
        self.typecButton.grid(row=2,column=0,sticky=tkinter.W,padx=35)
        self.typedButton = tkinter.Radiobutton(self.optionsframe, text="RRd", variable=self.RRLType, value=3)
        self.typedButton.grid(row=3,column=0,sticky=tkinter.W,padx=35)
        self.typeunButton = tkinter.Radiobutton(self.optionsframe, text="Unknown", variable=self.RRLType, value=4)
        self.typeunButton.grid(row=4,column=0,sticky=tkinter.W,padx=35)

        # Aperture check
        self.apOK = tkinter.IntVar()
        self.apOK.set(0)
        self.apOKLabel = tkinter.Label(self.optionsframe, text="Is Aperture OK?")
        self.apOKLabel.grid(row=0,column=1,sticky=tkinter.W,padx=30)
        self.apOKButton = tkinter.Radiobutton(self.optionsframe, text="OK", variable=self.apOK, value=1)
        self.apOKButton.grid(row=1,column=1,sticky=tkinter.W,padx=35)
        self.apMAYButton = tkinter.Radiobutton(self.optionsframe, text="Maybe", variable=self.apOK, value=2)
        self.apMAYButton.grid(row=2,column=1,sticky=tkinter.W,padx=35)
        self.apBADButton = tkinter.Radiobutton(self.optionsframe, text="BAD", variable=self.apOK, value=3)
        self.apBADButton.grid(row=3,column=1,sticky=tkinter.W,padx=35)

        # Period doubling region (yellow)
        self.peridoub = tkinter.IntVar()
        self.peridoub.set(0)
        self.peridoubLabel = tkinter.Label(self.optionsframe, text="Peak in yellow region?")
        self.peridoubLabel.grid(row=0,column=2,sticky=tkinter.W,padx=30)
        self.peridoubYESButton = tkinter.Radiobutton(self.optionsframe, text="YES", variable=self.peridoub, value=1)
        self.peridoubYESButton.grid(row=1,column=2,sticky=tkinter.W,padx=35)
        self.peridoubNOButton = tkinter.Radiobutton(self.optionsframe, text="NO", variable=self.peridoub, value=2)
        self.peridoubNOButton.grid(row=2,column=2,sticky=tkinter.W,padx=35)

        # Non-radial mode region (green)
        self.nonrad61 = tkinter.IntVar()
        self.nonrad61.set(0)
        self.nonrad61Label = tkinter.Label(self.optionsframe, text="Peak in green region?")
        self.nonrad61Label.grid(row=0,column=3,sticky=tkinter.W,padx=30)
        self.nonrad61YESButton = tkinter.Radiobutton(self.optionsframe, text="YES", variable=self.nonrad61, value=1)
        self.nonrad61YESButton.grid(row=1,column=3,sticky=tkinter.W,padx=35)
        self.nonrad61NOButton = tkinter.Radiobutton(self.optionsframe, text="NO", variable=self.nonrad61, value=2)
        self.nonrad61NOButton.grid(row=2,column=3,sticky=tkinter.W,padx=35)

        # Phased plot check
        self.phasedOK = tkinter.IntVar()
        self.phasedOK.set(0)
        self.phasedOKLabel = tkinter.Label(self.optionsframe, text="Is Phased Plot OK?")
        self.phasedOKLabel.grid(row=0,column=4,sticky=tkinter.W,padx=30)
        self.phasedOKButton = tkinter.Radiobutton(self.optionsframe, text="OK", variable=self.phasedOK, value=1)
        self.phasedOKButton.grid(row=1,column=4,sticky=tkinter.W,padx=35)
        self.phasedMAYButton = tkinter.Radiobutton(self.optionsframe, text="Maybe", variable=self.phasedOK, value=2)
        self.phasedMAYButton.grid(row=2,column=4,sticky=tkinter.W,padx=35)
        self.phasedBADButton = tkinter.Radiobutton(self.optionsframe, text="BAD", variable=self.phasedOK, value=3)
        self.phasedBADButton.grid(row=3,column=4,sticky=tkinter.W,padx=35)

        # Blazhko check
        self.isblazhko = tkinter.IntVar()
        self.isblazhko.set(0)
        self.isblazhkoLabel = tkinter.Label(self.optionsframe, text="Blazhko Effect?")
        self.isblazhkoLabel.grid(row=0,column=5,sticky=tkinter.W,padx=30)
        self.isblazhkoYESButton = tkinter.Radiobutton(self.optionsframe, text="YES", variable=self.isblazhko, value=1)
        self.isblazhkoYESButton.grid(row=1,column=5,sticky=tkinter.W,padx=35)
        self.isblazhkoMAYButton = tkinter.Radiobutton(self.optionsframe, text="Maybe", variable=self.isblazhko, value=2)
        self.isblazhkoMAYButton.grid(row=2,column=5,sticky=tkinter.W,padx=35)
        self.isblazhkoNOButton = tkinter.Radiobutton(self.optionsframe, text="NO", variable=self.isblazhko, value=3)
        self.isblazhkoNOButton.grid(row=3,column=5,sticky=tkinter.W,padx=35)

        # Aperture adjustment
        self.apselectLabel = tkinter.Label(self.apselectframe, text="Aperture Threshold:")
        self.apselectLabel.pack(side="left")
        self.apselectEntry = tkinter.Entry(self.apselectframe,width=4)
        self.apselectEntry.pack(side="left")
        self.apselectEntry.insert(0,'%.1f'%(self.apthresh))
        self.apupdateButton = tkinter.Button(self.apselectframe, text="Update", fg="green", command=lambda: self.showStar(True))
        self.apupdateButton.pack(side="left")

        # Additional notes
        self.notesLabel = tkinter.Label(self.notesframe, text="Notes:")
        self.notesLabel.pack(side="left")
        self.notesEntry = tkinter.Entry(self.notesframe,width=125)
        self.notesEntry.pack(side="left")

        # Aperture plot
        self.figContour, self.figAx = plt.subplots(figsize=(2,2))
        self.tpfPlot = FigureCanvasTkAgg(self.figContour, self.plotframe)
        self.tpfPlot.get_tk_widget().grid(row=0,column=0)

        # Lightcurve plot
        self.figLC, self.axLC = plt.subplots(figsize=(7,2))
        self.lcPlot = FigureCanvasTkAgg(self.figLC, self.plotframe)
        self.lcPlot.get_tk_widget().grid(row=0,column=1,columnspan=2)

        # Zoomed region of lightcurve
        self.figZoom, self.axZoom = plt.subplots(figsize=(3,2))
        self.zoomPlot = FigureCanvasTkAgg(self.figZoom, self.plotframe)
        self.zoomPlot.get_tk_widget().grid(row=0,column=3)

        # Lomb-Scargle analysis
        self.figPower, self.axPower = plt.subplots(figsize=(3,2))
        self.powerPlot = FigureCanvasTkAgg(self.figPower, self.plotframe)
        self.powerPlot.get_tk_widget().grid(row=1,column=0)

        # Pre-whitened Lomb-Scargle analysis
        self.figPreWhiten, self.axPreWhiten = plt.subplots(figsize=(3,2))
        self.prewhitenPlot = FigureCanvasTkAgg(self.figPreWhiten, self.plotframe)
        self.prewhitenPlot.get_tk_widget().grid(row=1,column=1)

        # Zoomed region of pre-whitened Lomb-Scargle analysis
        self.figPWZoom, self.axPWZoom = plt.subplots(figsize=(3,2))
        self.pwzoomPlot = FigureCanvasTkAgg(self.figPWZoom, self.plotframe)
        self.pwzoomPlot.get_tk_widget().grid(row=1,column=2)

        # Phased plot
        self.figPhased, self.axPhased = plt.subplots(figsize=(3,2))
        self.phasedPlot = FigureCanvasTkAgg(self.figPhased, self.plotframe)
        self.phasedPlot.get_tk_widget().grid(row=1,column=3)

        self.infoframe.pack(side="top",fill="x",pady=20)
        self.apselectframe.pack(side="top",fill="x")
        self.plotframe.pack(side="top",fill="x")
        self.optionsframe.pack(side="top",fill="x",pady=20)
        self.notesframe.pack(side="bottom",fill="x")

    ## Pre-whiten data by subtracting the
    ## primary frequency and next 9 harmonics
    def doPreWhiten(self):
        model_lc = self.lspg.model(self.star_lc.time*u.day,self.lspg.frequency_at_max_power)
        for i in range(2,10):
            model_lc = model_lc + self.lspg.model(self.star_lc.time*u.day,i*self.lspg.frequency_at_max_power) - 1.0
        diff2 = self.star_lc - model_lc + 1.0
        self.pwlspg = diff2.to_periodogram(maximum_frequency=12.0, oversample_factor=50)
        return

    ## Update number of stars left
    def changeCounter(self,N):
        self.fNumberLabel['text'] = "N Stars = %d"%(N)

    ## Update common name of star
    def changeTIC(self):
        self.sTICLabel['text'] = "%s"%(self.starName)

    ## Update RA,DEC
    def changeRADEC(self):
        self.sPosLabel['text'] = "(RA, DEC): (%f, %f)"%(self.RA[len(self.RA)-1],self.DE[len(self.DE)-1])

    ## Save options for star
    def storeValues(self):
        if(self.NDONE>=len(self.RTYPE)):
            return
        self.RTYPE[self.NDONE] = self.RRLType.get()
        self.APTEST[self.NDONE] = self.apOK.get()
        self.APSIZE[self.NDONE] = self.apthresh
        self.PHTEST[self.NDONE] = self.phasedOK.get()
        self.PDR[self.NDONE] = self.peridoub.get()
        self.NRR[self.NDONE] = self.nonrad61.get()
        self.BLA[self.NDONE] = self.isblazhko.get()
        self.ORA[self.NDONE] = self.RA[len(self.RA)-1]
        self.ODE[self.NDONE] = self.DE[len(self.DE)-1]
        self.NOTES[self.NDONE] = self.notesEntry.get()
        self.NAME[self.NDONE] = self.starName

        table = Table( [self.NAME, self.ORA, self.ODE, self.RTYPE, self.PERIOD, self.NSECS,
                        self.APTEST, self.APSIZE, self.PHTEST, self.PDR, self.NRR, self.BLA, self.NOTES],
                        names=('name','ra','dec','type','period','n sectors',
                                   'apcheck','apsize','phased','PD','NR','Blazhko','notes') )
        #table.write('tableRRL.fits', format='fits', overwrite=True)
        table.write(self.outfile, format='fits', overwrite=True)
        return

    ## Reset all options for star
    def clearValues(self):
        # reset information about the star
        return

    ## Popup a window to interactively look at the lightcurve
    def inspectLC(self):
        inspwin = tkinter.Tk()
        inspwin.title("Inspect Lightcurve")
        figInspLC, axInspLC = plt.subplots(figsize=(7,2))
        self.star_lc.scatter(ax=axInspLC)
        axInspLC.grid()
        lcInspPlot = FigureCanvasTkAgg(figInspLC, inspwin)
        #lcInspPlot.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(lcInspPlot, inspwin)
        toolbar.update()
        lcInspPlot.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

    def inspectPhased(self):
        phasedwin = tkinter.Tk()
        phasedwin.title("Inspect Phased Plot")
        t0 = self.star_lc.time[np.argmax(self.star_lc.flux)]
        figInspPhased, axInspPhased = plt.subplots(figsize=(6,3))
        self.star_lc.fold(period=self.lspg.period_at_max_power,t0=t0).scatter(ax=axInspPhased)
        yval = self.star_lc.flux[np.argmin(self.star_lc.flux)] + 0.2*(self.star_lc.flux[np.argmax(self.star_lc.flux)]-self.star_lc.flux[np.argmin(self.star_lc.flux)])
        axInspPhased.grid()
        phasedInspPlot = FigureCanvasTkAgg(figInspPhased, phasedwin)
        #phasedInspPlot.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(phasedInspPlot, phasedwin)
        toolbar.update()
        phasedInspPlot.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)


    ## Popup what is wrong
    ## DISABLED
    def errorWindow(self,errtext):
        return
#        errwin = tkinter.Tk()
#        errLabel = tkinter.Label(errwin, text=errtext)
#        errLabel.grid(row=0,sticky=tkinter.E)
#        errPosLabel = tkinter.Label(errwin, text="(RA, DEC): (%f, %f)"%(self.RA[len(self.RA)-1],self.DE[len(self.DE)-1]))
#        errPosLabel.grid(row=1,column=0)
#        errButton = tkinter.Button(errwin, text="OK", fg="red", command=lambda : errwin.destroy())
#        errButton.grid(row=2,column=0)

    ## Open an input file with target stars
    ## Format: ra(decimal degrees) dec(decimal degrees)
    def retrieveInput(self):
        fname = tkinter.filedialog.askopenfilename()
        endParse = fname.split("/")
        self.infile = endParse[-1]
        self.fNameLabel['text'] = "Input File Name: " + self.infile
        endParse = self.infile.split(".")
        self.outfile = endParse[0] + ".fits"
        temp = [0,1]
        self.RA.clear()
        self.DE.clear()
        try:
            for line in open(fname):
                parse = line.split()
                temp[0] = float(parse[0])
                temp[1] = float(parse[1])
                self.RA.append(temp[0])
                self.DE.append(temp[1])
            # Show number of stars in file
            self.changeCounter(len(self.RA))
            # We go backwards through the list
            self.RA.reverse()
            self.DE.reverse()
            self.PERIOD = np.zeros(len(self.RA), dtype='float')
            self.TESSID = np.zeros(len(self.RA), dtype='int')
            self.NSECS = np.zeros(len(self.RA), dtype='int')
            self.RTYPE = np.zeros(len(self.RA), dtype='int')
            self.APTEST = np.zeros(len(self.RA), dtype='int')
            self.APSIZE = np.zeros(len(self.RA), dtype='float')
            self.PHTEST = np.zeros(len(self.RA), dtype='int')
            self.PDR = np.zeros(len(self.RA), dtype='int')
            self.NRR = np.zeros(len(self.RA), dtype='int')
            self.BLA = np.zeros(len(self.RA), dtype='int')
            self.ORA = np.zeros(len(self.RA), dtype='float')
            self.ODE = np.zeros(len(self.RA), dtype='float')
            self.NOTES = np.empty(len(self.RA), dtype='U250')
            self.NAME = np.empty(len(self.RA), dtype='U100')
            self.NDONE = 0
            self.showStar()
        except:
            self.changeCounter(0)
            self.errorWindow("BAD INPUT FILE!!")

    ## Save options and move on
    def saveStar(self):
        self.storeValues()
        self.NDONE += 1
        if(len(self.RA)>0):
            self.RA.pop()
            self.DE.pop()
        self.RRLType.set(0)
        self.apOK.set(0)
        self.phasedOK.set(0)
        self.peridoub.set(0)
        self.nonrad61.set(0)
        self.isblazhko.set(0)
        self.notesEntry.delete(0, tkinter.END)
        self.changeCounter(len(self.RA))
        self.showStar()

    ## Plot the TPF image and aperture
    def addTPF(self):
        self.tpfPlot.get_tk_widget().grid_forget()
        #self.figContour.clear()
        plt.close(self.figContour)
        self.figContour, self.figAx = plt.subplots(figsize=(2,2))
        self.tpf[0].plot(ax=self.figAx, aperture_mask=self.target_mask0, mask_color='k')
        self.tpfPlot = FigureCanvasTkAgg(self.figContour, self.plotframe)
        self.tpfPlot.get_tk_widget().grid(row = 2, column = 0)

    ## Plot the full light curve and
    ## a zoomed in portion
    def addLC(self):
        self.lcPlot.get_tk_widget().grid_forget()
        #self.figLC.clear()
        plt.close(self.figLC)
        self.figLC, self.axLC = plt.subplots(figsize=(7,2))
        self.star_lc.scatter(ax=self.axLC)
        self.axLC.grid()
        self.lcPlot = FigureCanvasTkAgg(self.figLC, self.plotframe)
        self.lcPlot.get_tk_widget().grid(row = 2, column = 1, columnspan = 3)

        self.zoomPlot.get_tk_widget().grid_forget()
        #self.figZoom.clear()
        plt.close(self.figZoom)
        self.figZoom, self.axZoom = plt.subplots(figsize=(3,2))
        self.star_lc.scatter(ax=self.axZoom)
        self.axZoom.grid()
        self.axZoom.set_xlim(self.star_lc.time[0]+4,self.star_lc.time[0]+6)
        self.zoomPlot = FigureCanvasTkAgg(self.figZoom, self.plotframe)
        self.zoomPlot.get_tk_widget().grid(row = 2, column = 4)

    ## Plot a phased light curve
    def addPhased(self):
        # Get max value for t0 estimate
        t0 = self.star_lc.time[np.argmax(self.star_lc.flux)]
        self.phasedPlot.get_tk_widget().grid_forget()
        #self.figPhased.clear()
        plt.close(self.figPhased)
        self.figPhased, self.axPhased = plt.subplots(figsize=(3,2))
        # Fold the lightcurve using period from Lomb-Scargle analysis and t0
        self.star_lc.fold(period=self.lspg.period_at_max_power,t0=t0).scatter(ax=self.axPhased)
        # Print the period value on the plot
        yval = self.star_lc.flux[np.argmin(self.star_lc.flux)] + 0.2*(self.star_lc.flux[np.argmax(self.star_lc.flux)]-self.star_lc.flux[np.argmin(self.star_lc.flux)])
        self.axPhased.text(-0.2,yval,"P = %.5f"%(self.lspg.period_at_max_power.value))
        self.axPhased.grid()
        self.phasedPlot = FigureCanvasTkAgg(self.figPhased, self.plotframe)
        self.phasedPlot.get_tk_widget().grid(row = 3, column = 4)

    ## Plot analysis for an RRab/c type
    def addAnalysis_c(self):
        self.PERIOD[self.NDONE] = self.lspg.period_at_max_power.value
        fPrimary = self.lspg.frequency_at_max_power.value
        majorTicks = np.arange(0,13,1)
        minorTicks = np.arange(0,13,0.5)

        self.powerPlot.get_tk_widget().grid_forget()
        #self.figPower.clear()
        plt.close(self.figPower)
        self.figPower, self.axPower = plt.subplots(figsize=(3,2))
        # Plot Lomb-Scargle analysis
        self.lspg.plot(ax=self.axPower)
        # Print max frequency on plot
        self.axPower.text(1.1*fPrimary,0.8*self.lspg.max_power.value,"F = %.5f"%(fPrimary))
        self.axPower.set_xticks(majorTicks)
        self.axPower.set_xticks(minorTicks, minor=True)
        self.axPower.grid(which='minor', alpha=0.2)
        self.axPower.grid(which='major', alpha=0.5)
        self.powerPlot = FigureCanvasTkAgg(self.figPower, self.plotframe)
        self.powerPlot.get_tk_widget().grid(row = 3, column = 0)

        # Pre-whiten the Lomb-Scargle analysis
        self.doPreWhiten()

        self.prewhitenPlot.get_tk_widget().grid_forget()
        #self.figPreWhiten.clear()
        plt.close(self.figPreWhiten)
        self.figPreWhiten, self.axPreWhiten = plt.subplots(figsize=(3,2))
        # Plot pre-whitened Lomb-Scargle analysis
        self.pwlspg.plot(ax=self.axPreWhiten)
        self.axPreWhiten.set_xticks(majorTicks)
        self.axPreWhiten.set_xticks(minorTicks, minor=True)
        self.axPreWhiten.grid(which='minor', alpha=0.2)
        self.axPreWhiten.grid(which='major', alpha=0.5)
        self.prewhitenPlot = FigureCanvasTkAgg(self.figPreWhiten, self.plotframe)
        self.prewhitenPlot.get_tk_widget().grid(row = 3, column = 1)

        self.pwzoomPlot.get_tk_widget().grid_forget()
        #self.figPWZoom.clear()
        plt.close(self.figPWZoom)
        self.figPWZoom, self.axPWZoom = plt.subplots(figsize=(3,2))

        # Plot zoomed region of pre-whitened Lomb-Scargle analysis
        self.pwlspg.plot(ax=self.axPWZoom)
        # Print max frequency on plot
        self.axPWZoom.text(1.1*fPrimary,0.8*self.pwlspg.max_power.value,"F = %.5f"%(fPrimary))
        self.axPWZoom.set_xticks(majorTicks)
        self.axPWZoom.set_xticks(minorTicks, minor=True)
        self.axPWZoom.set_xlim(0,2.5*fPrimary)
        self.axPWZoom.grid(which='minor', alpha=0.2)
        self.axPWZoom.grid(which='major', alpha=0.5)
        # Gray regions of primary frequency and 1st harmonic
        self.axPWZoom.axvspan(0.97*fPrimary,1.03*fPrimary,alpha=0.42,color='gray')
        self.axPWZoom.axvspan(1.97*fPrimary,2.03*fPrimary,alpha=0.42,color='gray')
        # Green region for non-radial modes
        self.axPWZoom.axvspan(1.59*fPrimary,1.69*fPrimary,alpha=0.42,color='green')
        # Yellow region for period doubling
        self.axPWZoom.axvspan(1.45*fPrimary,1.55*fPrimary,alpha=0.42,color='yellow')
        # Red region for RRd mode
        self.axPWZoom.axvspan(0.72*fPrimary,0.78*fPrimary,alpha=0.42,color='red')
        self.pwzoomPlot = FigureCanvasTkAgg(self.figPWZoom, self.plotframe)
        self.pwzoomPlot.get_tk_widget().grid(row = 3, column = 2)

    ## Display all the plots
    def updatePlots(self):
        self.addTPF()
        self.addLC()
        self.addPhased()
        self.addAnalysis_c()

    ## Remove times with bad data
    ## NOT USED
    def goodTimes(self,lc,sec):
        if(sec == 20):
            return ((lc.time < 1856) | (lc.time > 1858))
        elif(sec == 17):
            return ((lc.time < 1772) | (lc.time > 1778))
        else:
            return (lc.time > 0)

    ## Get data and display
    def showStar(self,updateOnly=False):
        plt.close('all')
        ## Check to see if we've reached the end of the list
        if(len(self.RA)==0):
            self.changeCounter(0)
            self.errorWindow("ALL DONE.")
            self.changeCounter(len(self.RA))
            return

        self.apthresh = float(self.apselectEntry.get())
        
        ## Change RA/dec if this is a new star
        if(not updateOnly):
            ## clear marked values and update RA/DEC
            self.clearValues()
            self.changeRADEC()
            self.changeCounter(len(self.RA))

            ## get the coordinates of the next star
            coord = SkyCoord(self.RA[-1],self.DE[-1], unit="deg")
            ## get the name of this star from Simbad and update GUI
            try:
                self.starName = Simbad.query_region(coord,radius='0d0m2s')['MAIN_ID'][0]
            except:
                self.starName = "Unknown"
            self.changeTIC()

            ## download the TPFs for every possible sector
            self.tpf = search_tesscut(coord).download_all(cutout_size=(10,10))
            print(self.tpf)

        ## make sure the 1st sector is OK
        ntrgtpix = 0
        while(ntrgtpix == 0 and hasattr(self.tpf,"__len__") and len(self.tpf)>0):
            tmsk = self.tpf[0].create_threshold_mask(threshold=self.apthresh)
            ntrgtpix = tmsk.sum()
            if(ntrgtpix == 0):
                print('sector %d...SKIPPED.'%(self.tpf[0].sector), flush=True)
                self.tpf = self.tpf[1:]

        if(hasattr(self.tpf,"__len__") and len(self.tpf)>0):
            self.NSECS[self.NDONE] = len(self.tpf)
            print('sector %d...'%(self.tpf[0].sector), end='', flush=True)

            ## create an aperture mask for the target star
            self.target_mask0 = self.tpf[0].create_threshold_mask(threshold=self.apthresh)
            n_target_pixels = self.target_mask0.sum()

            ## create the lightcurve for the first sector using aperture
            raw_lc = self.tpf[0].to_lightcurve(aperture_mask=self.target_mask0)
            ## eliminate problem times from sector
            #time_mask = self.goodTimes(raw_lc,self.tpf[0].sector)
            #raw_lc = raw_lc[time_mask]

            ## create a mask for the background signal
            self.bkgr_mask0 = ~self.tpf[0].create_threshold_mask(threshold=0.001, reference_pixel=None)
            n_background_pixels = self.bkgr_mask0.sum()

            ## create a lightcurve for the background and scale it to the target aperture
            bkgr_lc = self.tpf[0].to_lightcurve(aperture_mask=self.bkgr_mask0) / n_background_pixels * n_target_pixels
            ## eliminate problem times from sector
            #bkgr_lc = bkgr_lc[time_mask]

            ## subtract the background from the signal
            self.star_lc = (raw_lc - bkgr_lc.flux)

            ## normalize the flux and remove problem points
            self.star_lc = self.star_lc.remove_nans().remove_outliers(sigma=6).normalize()
            #self.star_lc = self.star_lc[(self.star_lc.flux>0.6) & (self.star_lc.flux<1.6)]

            print('DONE.', flush=True)
            for itpf in self.tpf[1:]:
                print('sector %d...'%(itpf.sector), end='', flush=True)

                ## create an aperture mask for the target star
                target_mask = itpf.create_threshold_mask(threshold=self.apthresh)
                n_target_pixels = target_mask.sum()

                ## check to make sure the mask is OK
                if(n_target_pixels == 0):
                    print('SKIPPED.',flush=True)
                    continue

                ## create the lightcurve for the first sector using aperture
                raw_lc = itpf.to_lightcurve(aperture_mask=target_mask)
                ## eliminate problem times from sector
                #time_mask = self.goodTimes(raw_lc,itpf.sector)
                #raw_lc = raw_lc[time_mask]

                ## create a mask for the background signal
                bkgr_mask = ~itpf.create_threshold_mask(threshold=0.001, reference_pixel=None)
                n_background_pixels = bkgr_mask.sum()

                ## create a lightcurve for the background and scale it to the target aperture
                bkgr_lc = itpf.to_lightcurve(aperture_mask=bkgr_mask) / n_background_pixels * n_target_pixels
                ## eliminate problem times from sector
                #bkgr_lc = bkgr_lc[time_mask]

                ## subtract the background from the signal
                temp_lc = (raw_lc - bkgr_lc.flux)

                ## normalize the flux and remove problem points
                temp_lc = temp_lc.remove_nans().remove_outliers(sigma=6).normalize()
                #temp_lc = temp_lc[(temp_lc.flux>0.6) & (temp_lc.flux<1.6)]

                ## add this sector to the previous
                self.star_lc = self.star_lc.append(temp_lc)

                print('DONE.', flush=True)

            ## do a LombScargle analysis of the lightcurve
            self.lspg = self.star_lc.to_periodogram(maximum_frequency=12.0,oversample_factor=50)
        else:
            self.errorWindow("NOT IN SECTOR")
            self.saveStar()
            return

        ## Make plots
        self.updatePlots()
