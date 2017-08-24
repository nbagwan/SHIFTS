import os
import pdb
from scipy.stats.stats import linregress
import shutil
import numpy as np
import all_stats

separations = "/"

##### this scripts consists of two different method; 1) target_slope: its method to calculates gaussian derivatives withput using package
##### 2) SlopeCalculation, this method calculates a hitogram of all delta masss using intger binning by making making delta masses to center at 4th decimal point
# (to reduce the dispersion of masses###

#### the method starts with calculating the slope value between the masss bin and frequency with sliding window of 7 mass bins#####
#### next stpe is to calculate the slope value between the slope 1 and mass bins,  again with sliding window of 7 points #####

def target_slope(filelist, target_dic, bintoHist,target=True):
    outputlist = [
        ["Bin", "\t", "Frequency", "\t", "Slope1", "\t", "Slope2", "\t", "peak-Width", "\t", "peak-Apex", "\t",
         "intercept_mass", "\n"]]

    separations = "/"
    firstfilepath = os.path.dirname(filelist[0])
    if target:
        fileName=firstfilepath + "/target_Peak_identification_histogram.txt"
    else:
        fileName = firstfilepath + "/Decoy_Peak_identification_histogram.txt"

    w1 = open(fileName, "w")
    dictlist = []

    for key, value in target_dic.iteritems():
        temp = [key, value]
        dictlist.append(temp)

    sorted_dictlist = sorted(dictlist, key=lambda x: float(x[0]))
    binlist, frelist, looplist = [], [], []
    for p in sorted_dictlist:
        binlist.append(float(p[0]))
        frelist.append(float(p[1]))
        looplist.append(float(p[0]))

########## take the first 7 points and calculate the slope 1

    slope1 = [0]
    for index in range(0, len(looplist) - 6):
        s, intercept, r, p, std_error = linregress(binlist[index:index + 7], frelist[index:index + 7])
        slope1.append(s)

######### does the same thing but slope1 calculated above #####
    slope2 = []
    for index1 in range(0, len(looplist) - 12):
        s1, intercept1, r1, p1, std_error1 = linregress(binlist[index1 + 3:index1 + 10], slope1[index1 + 1:index1 + 8])
        slope2.append(s1)

    apex = []
    peak = []
    interceptList = [0]
    if len(looplist) % 2 == 0:
        minus1 = 6
        minus2 = 3
    else:
        minus1 = 7
        minus2 = 3

##### after slope 1 and 2 are calculated,  next loop checks for Apex by check he slope 1, (cehcks inflation) first value > 0 and next < 1

    for index3 in range(len(looplist) - minus1):
        if slope1[index3] > 0.0 and slope1[index3 + 1] < 0.0:
            apex.append("1")
        else:
            apex.append("0")

####### this loop tell the width of a slope by checking all values below 0
    for index4 in range(len(looplist) - 12):
        if slope2[index4] < 0:
            peak.append("1")
        else:
            peak.append("0")

    # print len(interceptList), len(slope1)
    slope1 = [0] * 2 + slope1 + [0] * 3
    slope2 = [0] * 6 + slope2 + [0] * 6
    apex = [0] * 3 + apex + [0] * (len(binlist) - (len(apex) + 3))
    peak = [0] * 6 + peak + [0] * 6

############ her in this loop program is calcuating a intercept between inflation and deflation point, which actually be the Peak apex.
    for index6 in range(len(looplist) - 6):
        if (abs(slope1[index6 + 1]) + abs(slope1[index6 + 2])) == 0.0:
            intercept_mass = float("inf")
            interceptList.append(intercept_mass)
        else:
            intercept_mass = looplist[index6] + (float(bintoHist) * abs(slope1[index6 + 1])) / (
                abs(slope1[index6 + 1]) + abs(slope1[index6 + 2]))
            #interceptList.append(all_stats.truncate(intercept_mass, 3))

   #### recently changed ######
            #interceptList.append(all_stats.check(intercept_mass, 4))
            interceptList.append(intercept_mass)
            #interceptList.append(float("%.3f" % intercept_mass) )

    interceptList = interceptList + [0] * (len(binlist) - len(interceptList))
    plot_x = []
    plot_y = []
    for index5 in range(len(looplist)):
        outputlist.append(
            [str(binlist[index5]), "\t", str(frelist[index5]), "\t", str(slope1[index5]), "\t", str(slope2[index5]),
             "\t", str(peak[index5]), "\t", str(apex[index5]), "\t", str(interceptList[index5]), "\n"])

        # print str(binlist[index5]), str(interceptList[index5])

    for write in outputlist:
        w1.writelines(write)
    w1.close()
    print "Finished with slope Calculations for target PSMs..... Copying files"


########## program write a output file,
    for eachfile in filelist[1:len(filelist)]:
        foldername = os.path.dirname(eachfile)
        shutil.copy(fileName, foldername)

    return fileName
##############################################################################################################################

##### SlopeCalculation method take the output of Comet_prioceesinScript and a bin to make histogram from user input.
##### example of bin = 0.001 ####

def SlopeCalculation(filelist, bintoHist):
    target_histDic = {}
    decoy_histDic = {}
    dictlist = []
    for file in filelist:
        with open(file) as PTMfile:
            next(PTMfile)
            for line3 in PTMfile:
                if line3 != "\n":
                    splits3 = line3.split("\t")
                    charge = int(splits3[2].strip())
                    corXcor = float(splits3[16].strip())
                    Delta_bin = float(bintoHist)
                    calibrated_Delta_MH = float(splits3[20].strip())
                    label = str(splits3[18].strip())

                    if calibrated_Delta_MH > 0:
                        intformula = int(calibrated_Delta_MH / Delta_bin) * Delta_bin + Delta_bin / 2 #### intformula makes the masses to center to 0.5 values using bins
                        #truncateDmass = float("%.3f" % intformula)
                        truncateDmass = all_stats.check(intformula, 4)
                        truncateDmass1 = truncateDmass

                    else:
                        intformula = int(calibrated_Delta_MH / Delta_bin) * Delta_bin - Delta_bin / 2
                        #truncateDmass = float("%.3f" % intformula)
                        truncateDmass = all_stats.check(intformula, 4)
                        truncateDmass1 = truncateDmass

                    if label == "Target":
                        if truncateDmass1 not in target_histDic:
                            target_histDic[truncateDmass1] = 1
                        else:
                            target_histDic[truncateDmass1] += 1
                    else:
                        if truncateDmass1 not in decoy_histDic:
                            decoy_histDic[truncateDmass1] = 1
                        else:
                            decoy_histDic[truncateDmass1] += 1


    #### finally two lists are created, one for target and one decoys and both are than send to target_slope method for guassian modelling ####
    targetFileName=["",""]
    targetFileName[0]=target_slope(filelist, target_dic=target_histDic, bintoHist=bintoHist,target=True)
    targetFileName[1]=target_slope(filelist, target_dic=decoy_histDic, bintoHist=bintoHist, target=False)
    return targetFileName