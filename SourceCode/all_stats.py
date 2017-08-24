import math
import numpy as np
import os

def master_table(masterFile):
    master_list = []
    child_list = []
    with open(masterFile) as master:
        next(master)
        for line in master:
            if line != "\n":
                splits = line.split("\t")
                if len(splits[0].strip()) != 0:
                    master_path = splits[0].strip()
                    if master_path not in master_list:
                        master_list.append(master_path)
                if len(splits[1].strip()) != 0:
                    child_path = splits[1].strip()
                    if child_path not in child_list:
                        child_list.append(child_path)
                else:
                    child_path = master_path
                    if child_path not in child_list:
                        child_list.append(child_path)

    return master_list, child_list


def correctedXcor(xcor, mass, charge):
    if charge < 3:
        R = 1.0
    else:
        R = 1.22
    CorXcor = math.log10(xcor / float(R)) / math.log10(2 * mass / 110)
    return CorXcor


def medianMethod(lst):
    if type(lst) == list:
        sortedLst = sorted(lst)
        lstLen = len(lst)
        index = (lstLen - 1) // 2

        if (lstLen % 2):
            return sortedLst[index]
        else:
            return (sortedLst[index] + sortedLst[index + 1]) / 2.0
    else:
        return lst


def truncate(f, n):
    return math.floor(f * 10 ** n) / 10 ** n

def check(mass, n):
    if mass < 0.0:
        massABS = abs(mass)
        testout = truncate(massABS, n)
        tempMass1 = str(testout).split(".")[1]
        tempMass0 = str(testout).split(".")[0]
        if len(tempMass1) < n:
            tempMass1new = tempMass1 + "0" * (n - len(tempMass1))
            output = "-" + tempMass0 + "." + tempMass1new
        else:
            output = "-"+str(testout)
    else:
        testout1 = truncate(mass, n)
        tempMass11 = str(testout1).split(".")[1]
        tempMass00 = str(testout1).split(".")[0]
        if len(tempMass11) < n:
            tempMass11new = tempMass11 + "0" * (n - len(tempMass11))
            output = tempMass00 + "." + tempMass11new
        else:
            output = str(testout1)
    return output

def createApexList(apexfile):
    apexmassList=[]
    for masline in open(apexfile):
        if masline != "\n":
            s = masline.split("\t")
            mass = float(s[0])
            apexmassList.append(float("%.6f" % mass))
            #apexmassList.append([truncate(mass, 3) for n in range(1)][0])
    return apexmassList


def getAlpha(alphafile):
    for line in open(alphafile):
        if line!="\n":
            if line.startswith("Alpha"):
                Alphavalue = line.split(":")[1].strip()
    return Alphavalue

def globalfdr(globalfdrfile):
    for line in open(globalfdrfile):
        if line!="\n":
            if line.startswith("globalFDR"):
                GlobalFDRvalue= line.split(":")[1].strip()
    return GlobalFDRvalue

def sigmaFinder(filename):
    for line in open(filename):
        if line!="\n":
            if line.startswith("threeSigma"):
                sigmavalue= line.split(":")[1].strip()
    return sigmavalue

def fullWidthHalfMaximum(MADfile):
    for line in open(MADfile):
        if line!="\n":
            if line.startswith("FWHM"):
                FWHM= line.split(":")[1].strip()
    return FWHM

def sigma2fwhm(sigma):
    return sigma * np.sqrt(8 * np.log(2))




def madCalculation(file, scoreInput):
    outputPath = os.path.dirname(file)
    madoutputFile = outputPath + "/MAD_and_FWHM_calculations.txt"
    writefile = open(madoutputFile, "w")
    madList = []
    abs_madlist = []
    SigmaFactor = 1.4826
    with open(file) as PTMfile:
        next(PTMfile)
        for line3 in PTMfile:
            if line3 != "\n":
                splits3 = line3.split("\t")

                label = splits3[18].strip()
                corXcor = float(splits3[16].strip())
                calibrated_delta_MZ = float(splits3[21].strip())
                xcorthershold = float(scoreInput)

                if label == "Target":
                    if corXcor >= xcorthershold and calibrated_delta_MZ > -0.1 and calibrated_delta_MZ < 0.1:
                        madList.append(calibrated_delta_MZ)

    sortedmadList = sorted(madList)

    for i in sortedmadList:
        abs_madlist.append(abs(i))

    median = np.median(abs_madlist)

    sigma = (SigmaFactor * median) ** 2

    MAD = math.sqrt(sigma)

    writefile.writelines("MAD:"+ "\t" + str(MAD) + "\n")
    # trunc_MAD = all_stats.truncate(MAD, 4)
    trunc_MAD = float("%.4f" % MAD)
    FWHM = sigma2fwhm(MAD)
    writefile.writelines("full width if peak at half maximum" + "\n")
    writefile.writelines("FWHM:" + "\t" + str(FWHM))
    writefile.close()

    return FWHM























