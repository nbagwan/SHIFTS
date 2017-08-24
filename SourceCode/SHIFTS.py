__author__ = "nbagwan"

"version = 20170823.0"

import glob
import os
import pdb
from optparse import OptionParser
import SlopeCalculation
import TaggingScript
import SlopeFDR
import SigmaCalculator
import PTMclassification
import Isotop_CorrectionAdvance
import Comet_PTM_processingScript
import PeakFDR
import all_stats
import time

start_time = time.time()

#usage = "usage: %mainMethod -P -B -X  -F -O"

usage = "usage:-P -B -X  -F -A -f -O -i\n Use -h for detailed help for parameter\nVersion = 20170823.0\nSHIFTS (Systematic, Hypothesis-free Identification of PTMs with controlled \nFDR " \
        "based on ultra-Tolerant database search)is a program made in the Jesus \nVazquez Cardiovascular Proteomics Lab at Centro Nacional de " \
        "Investigaciones \nCardiovasculares,for highthroughput PTM( Post translation modifications )processing."

parser = OptionParser(usage=usage)

parser.add_option("-P", "--Path2master", action="store", type="string", dest="Dir", help="Enter the path to masterFile")

parser.add_option("-B", "--BinSize", action="store", type="string", dest="binSize", help="Enter the bin size for peak modelling")

parser.add_option("-X", "--Xcor", action="store", type="string", dest="CorrXcor", help="Enter the Xcor thershold for Calibration")

parser.add_option("-F", "--Fastafile", action="store", type="string", dest="fastafile", help="Enter the fasta database")

parser.add_option("-A", "--ApexFilter", action="store", type="string", dest="apex", help="Enter an intiger for Frequncey filter")

parser.add_option("-f", "--FDRthreshold", action="store", type="string", dest="FDRs", help="Enter the FDR threshold; 0.01 for 1% FDR")

parser.add_option("-O", "--OutPutname", action="store", type="string", dest="Outname", help="Enter the output folder name")

parser.add_option("-i", "--isotopeCorrection", action="store", type="string", dest="isotop", help="Enter 1 if you want to use isotopcorrection else 0")


(options, args) = parser.parse_args()



if options.Dir:
    Dir = options.Dir
else:
    parser.error("Please enter the path to mater file -h")

if options.binSize:
    binSize = options.binSize
else:
    parser.error("Please enter the Bin size you want use -h")

if options.CorrXcor:
    CorrXcor = options.CorrXcor
else:
    parser.error("Please enter Xcor thershold you want use for Calibration -h")

if options.fastafile:
    fastafile = options.fastafile
else:
    parser.error("Please enter fasta file -h")

if options.apex:
    apex = options.apex
else:
    parser.error("Enter intiger for Frequncey filter")

if options.FDRs:
    FDRs = options.FDRs
else:
    parser.error("please enter FDR you want to use; 0.01 for 1% FDR")

if options.Outname:
    Outname = options.Outname
else:
    parser.error("Please enter name for output folder -h")

if options.isotop:
    isotop = options.isotop
else:
    parser.error("Please enter 0 for no, and 1 for yes")


separations = "/"

### "calibrate_every_path" is the main method where all the methods will be called as a forkflow for different calculation and processing ###
### under this method all other mwthod scripts will be called, such as "Comet_PTM_processingScript" for calibration. ###

def calibrate_every_path(pathTomasterFIle, bins, Xcorthershold, outputfolder, fasta, apexthersold, fdrFilter, isotopTF ):

    ##### the first method "Comet_PTM_processingScript.ProcessingFile" takes all the raw comet output from master file and calculates a machine error and creates an output folder###
    #### for more details check the method Comet_PTM_processingScript.py#####

    filelist = []
    master_list, child_list = all_stats.master_table(masterFile=pathTomasterFIle)
    for subfolder_path in child_list:
            print "processing ", subfolder_path
            filelist += [Comet_PTM_processingScript.ProcessingFile(folder=subfolder_path, outfoldername=outputfolder, scoreInput=Xcorthershold)]


#####################################################################################################################################################################
########## this method calculates the delta peak width for PSM to peak apex assigment #####
    print "Starting with three sigma calculations for delta peak width for PSM to peak apex assignment"
    SigmaCalculator.filemakerforSigmaFinder(scoreInput=Xcorthershold, processingFileList=filelist, fdrthershol=fdrFilter)

#####################################################################################################################################################################
####### this method calulates a calculates the PEak paexes from a histogram of delta masses #####
    ##### for more details check the program code ####
    print "started with slope calculation, sliding window for slope calculation is 7."
    histogramFilenameList = SlopeCalculation.SlopeCalculation(filelist=filelist, bintoHist=bins)
    print "Slope calculation part DONE"
    print("---%s seconds ---" % (time.time() - start_time))

#####################################################################################################################################################################
########## this is a method for local FDR Calculations ################
    print "starting the slope FDR calculation......"
    SlopeFDR.slopeFDR(processingFileList=filelist)
    print("---%s seconds ---" % (time.time() - start_time))

#####################################################################################################################################################################
######## this a method for creat a list of peak apex and tag every PSMs to peak apex ######
    pathListforPTMclassifier = []
    for calibrationfile in filelist:
        filepath = os.path.dirname(calibrationfile)
        pathListforPTMclassifier.append(filepath)

        apexfile_T = TaggingScript.MassTagging2Seq_T(histogramFilenameList[0], apexthersold, target=True)
        if histogramFilenameList[1] == "":
            apexfile_D = ""
        else:
            apexfile_D = TaggingScript.MassTagging2Seq_T(histogramFilenameList[1], apexthersold, target=False)

#####################################################################################################################################################################

        subfolder = os.path.dirname(filepath)
        fastaFile = fasta
        slopeFilename = filepath + "/" + "SlopeFDRfile.txt"
        sigmaFile = filepath + "/" + "sigmaCalculations.txt"
        globalFDRfile = filepath + "/" + "globalFDRforSmallDB.txt"

        all_stats.madCalculation(file=slopeFilename, scoreInput=Xcorthershold)

        madfiles = filepath + "/" + "MAD_and_FWHM_calculations.txt"

        PeakFDR.peak_fdr_method(SlopeFDRfile=slopeFilename, apexList= apexfile_T, FWHM= madfiles)
        peakfile = filepath + "/Peak_and_Slope_FDRfile.txt"



        TaggingScript.taggingMethodForCOMET_PTM([apexfile_T, apexfile_D], calibrationfile, filepath,
                                                fasta=fastaFile, slopefdrFile=peakfile, globalfdrfile=globalFDRfile, sigmafile=madfiles, fdrTOuse=FDRs)



##### this a method to classify PSMs which can not be tagged to in method abouve ##################
    PTMclassification.ptmClassifier(listofpaths=pathListforPTMclassifier, sigmaFIle="sigmaCalculations.txt")

    if str(isotopTF) == "1":
        for everyPath in pathListforPTMclassifier:
            isotopFIle = Isotop_CorrectionAdvance.isotop(targetFile="AllWithSequence-massTag.txt", folder=everyPath, madFIle=madfiles)


if __name__ == "__main__":
    calibrate_every_path(pathTomasterFIle=Dir, Xcorthershold=CorrXcor, outputfolder=Outname, bins=binSize, fasta=fastafile, apexthersold=apex, fdrFilter=FDRs, isotopTF=isotop)

print("---%s seconds ---" % (time.time() - start_time))
    #
    # PDver="2.1"
    # folder="E:/Users/nbagwan/Desktop/isPTM_FDR_Comet/Test_FDR/masterFile_PESA_Experiment.txt"
    # outfoldername="PESA_FDR_test_27"
    #
    #
    # calibrate_every_path(PDversion=PDver, SmallDatabasePresent="0",pathTomasterFIle=folder, Xcorthershold="0.15", Experiment="TMT10",
    #                          outputfolder=outfoldername, quantificationPresent="0", bins="0.001", slidingWindow= "7",
    #                      fasta="Human_jul14_concate.fasta", apexthersold="0")
