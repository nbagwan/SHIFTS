# **Systematic, Hypothesis-free Identification of PTMs with controlled FDR based on ultra-Tolerant database search (SHIFTS)**
# **Manual, build (20170823.0)**

### Introduction

**SHIFTS** is a program made in the Jesus Vazquez Cardiovascular Proteomics Lab at Centro Nacional de Investigaciones Cardiovasculares (CNIC), for high throughput PTM (Post translation modifications) processing
SHIFTS identifies peaks in the Delta Mass distribution, assigns PSM to peaks and calculates FDR for peptide identification. 
SHIFTS is developed in python programming language as command line tool and its works with cometPTM produced “.txt” files.

### Availability

[SHIFTS is available at](https://github.com/nbagwan/SHIFTS.git)
Download the program folder in local to use, in program folder there are two sub folders; Dist and Build. The executable for SHIFTS is available in Dist folder. From here tool can be accessed, and for detailed help for parameter type

SHIFTS.exe -h

### Parameters

Usage: -P -B -X -F -A -f -O –i
Use –h or --help for detailed help for parameter

There are 8 input parameters for SHIFTS

commond line input	| explaintions
--------------------| ------------
-P, --Path2master	| Input will be a text file, containing two columns: first would be master path and second as subpath.this parameter will help to analyze multiple experiment of a project together, i.e Project name = PTM , experiments in project PTM; heart, liver. So, path/PTM will be masterPath and heart and liver will be subpath
-B, --BinSize	| Bin size for slope modelling, the most standard bin to use is 0.001 which has been tested for several experiments, however you are free to choose depending on your machine resolution. Lower the resolution bigger the bin
-X, --Xcor		| Corrected Xcorr threshold for choosing the best non-modified peptides, based on provided Xcorr, mass calibration will be performed. In most cases 0.20 or 0.25 would be good enough
-F, --Fastafile	| path to Protein database which was used for searches (Concatenated)
-A, --ApexFilter	| Apex threshold, a simple integer input (example 10 0r 15 0r 8) for ignoring the PTM peaks which have less PSMs than specified number. It will help to avoid the background noise. However in cases of small experiments using 0 would be useful, as some good PTMs will be very low frequent
-f, --FDRthreshold	| FDR filtration threshold. Example: 0.01 0r 0.05
-O, --OutPutname	| output folder name, specified name will be created in every experiment folder with all the results file
-i, --isotopeCorrection	| it´s True/false condition for isotopic correction. 0 for false and 1 for true. A simple isotopic correction to minimize misassignations of the correct monoisotopic peak of the precursor. When two PSM having the same sequence are encountered having a ?Mass difference within 1 ppm of the mass difference expected for either one or two 13C or one 34S, the ?Mass of the heaviest precursor is substituted by that of the lighest one


example to commond line. 
**SHIFTS.exe –P Path:\masterTable.txt -B 0.001 -X 0.22 -F Path:\gencode25.TargetDcoyConcat.fasta  -A 0 -f 0.01 -O SHIFTS-output -i 0**

### Output files
file name	| explaintions
------------| ------------
"experimentName_TargetData_Calibration.txt" | This is the first file to be created by SHIFTS, will contain all the fraction from a experiment with the calibration error per fraction with recalibrated masses. from here to the final results this file will be used
"GlobalFDR.txt"| This file will contain the global FDR Xcorr threshold for specified FDR in input parameter "-f"
"sigmaCalculations.txt"| This file will contain the fitting value for sigma
"heart_Median_log.txt"| This file contain all the calibration error value, and number of PSM/Fraction from which error has been calculated
"target_Peak_identification_histogram.txt"| This file contains a histogram and guassian derivatives for all the delta masses, which are used for peak apex picking
"Decoy_Peak_identification_histogram.txt"| This file contains a histogram fro decoys
"SlopeFDRfile.txt"|This file contains all the data, like the calibration file, but with Local FDR calculated with every 1 Da
"Peak_and_Slope_FDRfile.txt"|This file contains all the data, like the SlopeFDr file, but with Peak FDR calculated with every peak apex
"NotassignedSequences.txt"| This file contains some duplicated (for computing convenience) information whcih will be found in final file
"AllWithSequence-massTag.txt"| This file contains all the data with FDR filtration and PSM to peak apex tagging.(PSM to Peak assignement). this will be a final File in case of where -i from parater defined as 0
"IsotopCorrection_TargetData_withSequence-massTag.txt"| Will be the final file in case of where -i is defined as 1. contains all the PSMs after FDR filter and with Isopic corrections