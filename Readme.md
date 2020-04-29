# TRLFS analysis tool
This tool can be used to plot x,y data, perform simple background subtraction and integrate a data series.  
It comes with a GUI to facilitate usage.
## Basic usage
Input can be either a text file such as .asc, .txt, .dat or the Andor Software File .sif. Multi-imaged .sif files,
especially lifetime data will be transferred into individual files from one .sif file on opening.

###Open    
Opens a tkinter dialog for opening one or multiple files, which is/are then added to the lister box below.
Data of all files is plotted.

*Features*
* Files can be removed from the box by selecting them and pressing the 'del' button on the keyboard.
* .sif files can be read and, if multi-image .sif files are used it will divide them into individual .asc files which
can be saved as such
* remembers last folder from which files were loaded for quicker navigation

###Baseline subtraction    
Opens a tkinter window with the files selected in the lister box.
Double-click adds draggable nodes for the baseline creation, right-mouse button removes a selected point.

###Choose integration range    
Gives two draggable lines to define an integration range.

###Save integration range
Saves x,y data if data was integrated before by using either **Integreate all** or **Integrate**

###Save selected data
Saves x,y data as .asc file for the data selected in the lister box in a selected folder.

###Save all
Saves all files in listerbox as .asc files in selected folder.

###Plot
Plots selected data.    
*Note*: Double-clicking an item from the lister box or selecting multiple files using 'Control' or 'Shift' and pressing 
'Enter' will do the same.

###Clear list
Removes all files

###Integrate   
Integrates the selected files and presents the resulting values.

*Feauters*  
*   If files are derived via opening a multi-image .sif file, the delay step of the individual files is the x-axis.
*   If files end on a wavelenght in the shape of abc_58000.asc (in case of a wavelength of 580.00 nm) the wavelenght
will be uesd as x-axis

###Integrate all
Integreates all files in the lister box and presents the resulting values.

*Feauters*  
*   If files are derived via opening a multi-image .sif file, the delay step of the individual files is the x-axis.
*   If files end on a wavelenght in the shape of abc_58000.asc (in case of a wavelength of 580.00 nm) the wavelenght
will be uesd as x-axis

## Getting started
Besides input x,y data (or a .sif file) of any kind you will only need a few libraries:
* numpy
* scipy
* matplotlib
* tkinter
* trlfs (included in this repository)
* sif_reader (fujiisoup)

The sif_reader has to be modified to allow full support of .sif Andor files by overwriting the _sif_open.py file in
/venv/Lib/site-packages/sif_reader/ with the updated file in this repository.

If all is set just run the TRLFS_analysis_tool.py to start up the GUI.

## Try it out
The repository comes with a folder called examples, which includes a multi-image .sif file such as some .asc files.
When all libraries are 

## License statements
The license of sif_reader includig the _sif_open.py file redistributed here underlies the following terms:

Copyright (c) 2006, Marcel Leutenegger All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the 
following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
 in the documentation and/or other materials provided with the distribution
* Neither the name of the Ecole Polytechnique Fédérale de Lausanne, Laboratoire d'Optique Biomédicale nor the names of its
 contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

