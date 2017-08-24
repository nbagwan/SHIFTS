__author__ = "nbagwan"
import pdb
import os
import all_stats

####### this is method to calculate peak FDR for every peak apex, using the sigma value #######


def peak_fdr_method(SlopeFDRfile, apexList, FWHM):
    takeClosest = lambda num, collection: min(collection, key=lambda x: abs(x - num))
    # for calibrationfile in calibrationFilepaths:
    mad_fwhm= all_stats.fullWidthHalfMaximum(FWHM)
    filepath = os.path.dirname(SlopeFDRfile)
    #     slopefdrFile = filepath + "/" + "SlopeFDRfile.txt"
    outputfilename = filepath + "/Peak_and_Slope_FDRfile.txt"
    w = open(outputfilename, "w")
    fnn_sc_chDic = {}
    Slope_peak_Dict = {}
    outputFilelist = ["Scan\tSearchengineRank\tcharge\texpMH\ttheoMH\tExpmz\tXcor\tSeq\tRetentionTime\tProtAcc\tDeltaMod"
                 "\tB-series\tY.series\tIsotpicJump\tdeltaPeptide\tfilename\tCorXcor\tnewExpMH\tLabel\tMedian"
                 "\tCal_Delta_MH\tCal_Delta_M/Z\tCalExp_MZ\t1DawindowCenter\tdecoyCount\ttargetCount\tlocalFDR\t1Da_start"
                      "\t1Da_end\tPeakApex\tpeakDecoyCount\tpeakTargetCount\tPeakFDR\n"]

    outputFilelist_dic = {}
    with open(SlopeFDRfile) as slopeFile:
        next(slopeFile)
        for ln in slopeFile:
    #for ln in open(SlopeFDRfile):
            if ln != "\n":
                newsplits = ln.split("\t")
                sc = newsplits[0].strip()
                ch = newsplits[2].strip()
                xc = newsplits[16].strip()
                fnn = newsplits[15].strip().split(".")[0]
                fdr = newsplits[26].strip()
                st = newsplits[27].strip()
                ed = newsplits[28].strip()
                slopePeak = newsplits[23].strip()
                st_end = st + "_" + ed
                fnn_sc_ch = fnn + ".raw" + "_" + sc + "_" + ch
                if float(slopePeak) not in Slope_peak_Dict:
                    Slope_peak_Dict[float(slopePeak)] = [ln.strip()]
                else:
                    Slope_peak_Dict[float(slopePeak)].append(ln.strip())


    slopePeak_list = Slope_peak_Dict.keys()
    ApexList = []
    for apexline in open(apexList):
        if apexline not in ApexList:
            ApexList.append(float(apexline.strip("\n")))

    for everyApex in ApexList:
        closeinSlopePeak_List = takeClosest(everyApex, slopePeak_list)
        if closeinSlopePeak_List in Slope_peak_Dict:
            for element in Slope_peak_Dict[closeinSlopePeak_List]:
                theo_mass = float(everyApex)
                experimental_mass = float(element.split("\t")[20])
                corrXcor = element.split("\t")[16]
                label = element.split("\t")[18]
                charge = element.split("\t")[2]

                ppm_error = abs((experimental_mass - theo_mass ) * 1000000 / theo_mass)
                massDiff = abs(experimental_mass- theo_mass)
                #if ppm_error <= 10:
                if massDiff <= float(charge) * float(mad_fwhm) : ##### checking if a scan falls within a PTM peak by using MAd and Charge
                #if massDiff <= float(mad_fwhm):
                    #outputFilelist.append(element)
                    if everyApex not in outputFilelist_dic:
                        outputFilelist_dic[everyApex] = [[float(corrXcor), label, element, str(everyApex)]]
                    else:
                        outputFilelist_dic[everyApex].append([float(corrXcor), label, element, str(everyApex)])

                    Slope_peak_Dict[closeinSlopePeak_List].remove(element)


    for apex_mass in outputFilelist_dic:
        outputFilelist_dic[apex_mass].sort(key=lambda row: (row[0]), reverse=True)
        decoy, target = 0.0, 0.0
        for everyline in outputFilelist_dic[apex_mass]:
            if everyline[1] == "Decoy":
                decoy += 1
                everyline.append(decoy)
                everyline.append(target)
                try:
                    PeakFDR = float(decoy / target)
                except ZeroDivisionError:
                    PeakFDR = "1"
                everyline.append(PeakFDR)
                outputFilelist.append("\t".join(everyline[2].split("\t")) + "\t" + str(everyline[3]) + "\t"+ str(decoy) + "\t"+ str(target) + "\t" + str(PeakFDR) + "\n" )

            elif everyline[1] == "Target":
                target += 1
                everyline.append(decoy)
                everyline.append(target)
                try:
                    PeakFDR = float(decoy / target)
                except ZeroDivisionError:
                    PeakFDR = "1"
                everyline.append(PeakFDR)
                outputFilelist.append("\t".join(everyline[2].split("\t")) + "\t" + str(everyline[3]) + "\t"+ str(decoy) + "\t" + str(target) + "\t" + str(PeakFDR) + "\n")

    for restMass in Slope_peak_Dict:
        for mass in Slope_peak_Dict[restMass]:
            outputFilelist.append("\t".join(mass.split("\t")) + "\t" + "NA" + "\t" + "1"+ "\t" + "1"+ "\t" + "1" + "\n")

    for towrite in outputFilelist:
        w.write(str(towrite))
    w.close()
