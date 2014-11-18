------------------------------------
Analysing the data in Gaspware-XTRACKN
------------------------------------

For Beginners:

	Click R 	- open an energy spectra in the format: G0#D.000X|l:16        
	                   X is the run number and D the detector number. Start with D = 00 
	       View the different detectors: left click: [#+] or [#-] buttons in the low left corner 
	       View the different runs:     right click: [#+] or [#-] buttons in the low left corner
	
	
	Click Open CM   - open a compressed matrix .CMAT (the Y projection)
	                 
               W  =  set gates
               CW =  display the gated projection 
               ZW =  delete gates 
               MW =  display gates
               DW =  list/delete/group gates 
               Q  =  display initial projection of the matrix

	
	
	Recalibrate spectra (0.5kev/chan) click [EnCal] -> n (new) -> A(1)=0.5 -> (leave the rest empty)
	
	Zoom in X: put the markers with SPACE and type 'e'(expand)
	Zoom out: type 'ff'
	Zoom in Y: mousewheel
	Show peak energies: 'cp'
	Find a specific energy, type 'p' and enter the energy value.
	
	
	Auto-fit individual peaks: hold CTRL and left click the peak; in the terminal the energy, area and FWHM will be displayed
	Fit (integrate) multiple peaks: set fit region with 'r', set peaks with 'g', select background left+rigt with 'b'. Click 'Fit'
	Integrate manually: set integration region with 'i', background region with 'b' and click INT (middle right part of xtrackn)
	!!! Before each fitting procedure type 'za' to delete all the previous markers
	
	
        Clear the screen: type '='.
	To add aditional spectra: Right Click and select "Add Line" or "Add Column"
	To cancel any command added in the terminal type "Ctrl+D"
	To exit: "Ctrl+C"
	
DF/ZF = open/close file for storing fit data
DK    = Energy calibration / Peak-search signifiance 
CV    = Fit Gaussian
CJ    = Integrate
SX/SY = Same X / Same Y



      
	  How to calibrate the energy with XTrackn


 For HPGe

1. Load calibration spectra
2. Fit (ctrl+click) 2 peaks and click Cal2P. Input the correct energies
      If 'cp' does not identify all peaks, modify peaksearch significance with 'dk' + enter, enter, ... or select EnCal and input by channel
3. Click DT. Select calibration source. Enter, Enter... Select polynomial of 3rd/4th degree for a better fit. 
4. Save coeff in file. 
5. Move to the next detector (#+)
6. Ctrl+Click(DT) - will perform automatically step 4
7. Repeat to infinity from step 5.


FOR LaBR: at step 4 also input FWHM fit region (fit a peak after Cal2P to find the approx FWHM)

Example of custom file containing calibration peaks:

calibrate_planar_GE.txt with X-rays from 152Eu: (add energies by hand or by file)
           
             39.50000	   2.9999999E-02
             40.11000	   2.9999999E-02
             45.30000	   2.9999999E-02
             46.57000	   2.9999999E-02
             121.7817	   3.0000001E-04
             244.6975	   7.9999998E-04
             344.2785	   1.7000000E-03
             411.1165	   8.0000004E-03
             443.9650	   6.0000001E-03
     
           
           
        "Run by Run" calibration:
        
1. First calibrate with Eu (see above steps)
2. Create file with gammas from the runs you want to calibrate (file should have 2 columns: Energy Error )
3. DT -> select from file
4. After you calibrate de first detector, type CT
5. For filename input =.+ (mantains the same name but increments the extension)
6. Many enters for default options.
7. Wait command (...) [N] -> type Y. Otherwise it will do automatically everything and you will skip the magic.  
8. At the end of the runs, CTRL+D and change to the next detector. Repeat from step 4.
9. The recal command for gsort should contain RUN specification:
	  RECAL G0 Ge_ener.cal RUN 0 2 10 2047 2 20

	  
	  
