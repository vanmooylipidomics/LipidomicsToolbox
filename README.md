# LipidomicsToolbox (TESTING)
R scripts for converting &amp; processing HPLC-MS lipid data. These scripts were developed in the [Van Mooy Lab](http://www.whoi.edu/page.do?pid=80356) at [Woods Hole Oceanographic Institution](http://www.whoi.edu/). Many of them were written specifically to prepare HPLC-ESI-MS data from an Exactive Orbitrap mass spectrometer for follow-on analysis with the [LOBSTAHS](http://github.com/vanmooylipidomics/LOBSTAHS) package.

<h4>New R scripts in this repository:</h4>
   
 1. [prepOrbidata_shinydemo.R](https://github.com/vanmooylipidomics/LipidomicsToolbox/blob/testing/prepOrbidata_shinydemo.R):  TESTING. This is nearly identical to prepOrbidata.R but has been refactored for XCMS3 and includes a GUI option for entering peak-picking, RT correction, and grouping steps. Intended to use as a calibration tool for specific MS experiments.

 2. [StandardCheckApp.R](https://github.com/vanmooylipidomics/LipidomicsToolbox/blob/testing/StandardCheckApp.R):  TESTING. Tool for calibrating MS stadard compounds in a GUI. Will screen for expected compounds (standards) at a given m/z range and return abundance, RT, and can plot chromatogram of the found peaks in the experiment. Must run peak-picking first before using.
 
 
 
<h4>Data (.mzXML):</h4>

The raw data files that used to reside here have been incorporated into their own R package **PtH2O2lipids** and thus have their own repository: https://github.com/vanmooylipidomics/PtH2O2lipids

<h4>A note on CAMERA and xcms versus other tools</h4>

While the scripts here use [xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html) and [CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html) to accomplish much of the necessary processing (peak picking, retention time correction, grouping, etc.), we've also used the fast and very good program [MAVEN](http://genomics-pubs.princeton.edu/mzroll/index.php) for many projects. In fact, the older version of our analysis pipeline relied on MAVEN for peak picking, chromatographic alignment, grouping, and assignment of compound identities from our databases; the scripts for the pipeline are still in the repository https://github.com/vanmooylipidomics/old_pipeline.
