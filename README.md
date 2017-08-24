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

###Parameters
Usage: -P -B -X -F -A -f -O –i
Use –h or --help for detailed help for parameter

There are 8 input parameters for SHIFTS

commond line input	| explaintions
--------------------| ------------
-P, --Path2master   | Input will be a text file, containing two columns: first would be master path and second as subpath.this parameter will help to analyze multiple experiment of a project together, i.e Project name = PTM , experiments in project PTM; heart, liver. So, path/PTM will be masterPath and heart and liver will be subpath
					