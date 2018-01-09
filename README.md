# **Systematic, Hypothesis-free Identification of PTMs with controlled FDR based on ultra-Tolerant database search (SHIFTS)**
# **Manual, build (20170823.0)**

### Introduction

**SHIFTS** is a program made in the Jesus Vazquez Cardiovascular Proteomics Lab at Centro Nacional de Investigaciones Cardiovasculares (CNIC), for high throughput PTM (Post translation modifications) processing
SHIFTS identifies peaks in the Delta Mass distribution, assigns PSM to peaks and calculates FDR for peptide identification. 
SHIFTS is developed in python programming language as command line tool and its works with cometPTM produced “.txt” files.

### Availability

[SHIFTS is available at](SHIFTS_executable)
Download the program folder in local to use, in program folder there are two sub folders; Dist and Build. The executable for SHIFTS is available in Dist folder. From here tool can be accessed, and for detailed help for parameter type

SHIFTS.exe -h

### Parameters

Usage: -P -B -X -F -A -f -O –i -l -p
Use –h or --help for detailed help for parameter

There are 8 input parameters for SHIFTS

commond line input	| description
--------------------| ------------
-P, --Path2master	| Input will be a text file, containing two columns: first would be master path and second as subpath.this parameter will help to analyze multiple experiment of a project together, i.e Project name = PTM , experiments in project PTM; heart, liver. So, path/PTM will be masterPath and heart and liver will be subpath
-B, --BinSize	| Bin size for slope modelling, the most standard bin to use is 0.001 which has been tested for several experiments, however you are free to choose depending on your machine resolution. Lower the resolution bigger the bin
-X, --Xcor		| Corrected Xcorr threshold for choosing the best non-modified peptides, based on provided Xcorr, mass calibration will be performed. In most cases 0.20 or 0.25 would be good enough
-F, --Fastafile	| path to Protein database which was used for searches (Concatenated)
-A, --ApexFilter	| Apex threshold, a simple integer input (example 10 0r 15 0r 8) for ignoring the PTM peaks which have less PSMs than specified number. It will help to avoid the background noise. However in cases of small experiments using 0 would be useful, as some good PTMs will be very low frequent
-f, --FDRthreshold	| Global FDR filtration threshold. Example: 0.01 0r 0.05
-O, --OutPutname	| output folder name, specified name will be created in every experiment folder with all the results file
-i, --isotopeCorrection	| it´s True/false condition for isotopic correction. 0 for false and 1 for true. A simple isotopic correction to minimize misassignations of the correct monoisotopic peak of the precursor. When two PSM having the same sequence are encountered having a ?Mass difference within 1 ppm of the mass difference expected for either one or two 13C or one 34S, the ?Mass of the heaviest precursor is substituted by that of the lighest one
-l, --localFDR		| Local FDR filtration threshold. Example: 0.01 0r 0.05
-p, --localFDR		| Peak FDR filtration threshold. Example: 0.01 0r 0.05

example to commond line.

**SHIFTS.exe –P Path:\masterTable.txt -B 0.001 -X 0.22 -F Path:\gencode25.TargetDcoyConcat.fasta  -A 0 -f 0.01 -l 0.01 -p 0.01 -O SHIFTS-output -i 0**


### Output files
file name	| description
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

### Information and columns in output file

Thr final output file contains 39 columns including all possible information from every step in SHIFTS, belonging to every scan
Some important column which were added on top of cometPTM output in SHIFTS are

1. If -i in input parameter is set to 0
	1. CorXcor: This is score SHIFTS uses, which is based on Xcorr from cometPTM. takes in account the charge and average mass of amino acid. Input parameter -X uses this score
	1. median: This is median value in m/Z which reflects the machine error and is used for mass calibration
	1. "Calibrated DeltaMH",	"Calibrated Delta MZ",	"Calibrated EXP MZ": contains calibrated values
	1. SlopeFDR: contains the local FDR value
	1. TagsforPTM-Orphans: tag for PSMs as Orphan or NA for PTM or Non orphan Peak apex
	1. FinalSeq_Mass: contains the sequence provided with cometPTM, with one change that mass in sequence is now peak apex

1. If -i in input parameter is set to 1, these column will be additional
	1. Corr_Seq-mass: sequence with peak apex, after isotopic correction
	1. Corr_mass: just the delta mass, after isotopic correction
	1. Monoisotop_T/F: this column is indicated if isotopic correction was done for a scan(Yes/No)
	
### Example data set
[EXAMPLE can be foud here](EXAMPLE)
	
1. Download folder "exampleData.rar" and extract, which contains the two experiment "Heart" and "Liver". 
1. Download the MasterFile.txt as well and change "path" in masterPath column path/Project-mito and subpath column path/Heart.
1. Download the SHIFTS executable folder and on command line, change directory to "dist"
1. Run program as explained above.

