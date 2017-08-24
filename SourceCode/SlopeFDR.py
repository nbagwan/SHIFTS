__author__ = "nbagwan"
import pdb
from scipy.stats.stats import linregress
import os
import shutil
import all_stats

### this is a method for calculating the local fdr for evry approx. 1 Da window. its takes the input from the first step Comet_processingScript.


def slopeFDR(processingFileList):

    #alpha = all_stats.getAlpha(alphafile)
    alpha = 1
    for file in processingFileList:
        print "Started with slopeFDR for:", os.path.dirname(file)
        with open(file) as PTMfile:
            slopeMassList = []
            slope_intercept = []
            slopeMedian = []
            next(PTMfile)
            for line3 in PTMfile:
                if line3 != "\n":
                    splits3 = line3.split("\t")
                    calibrated_Delta_MH = float(splits3[20].strip())
                    slopeMassList.append(calibrated_Delta_MH)
                    if calibrated_Delta_MH > 0:
                        intOFmass = int(calibrated_Delta_MH + 0.5)
                    else:
                        intOFmass = int(calibrated_Delta_MH - 0.5)

                    slope_intercept.append(intOFmass)

            minOFslope = min(slope_intercept)
            maxOFslope = max(slope_intercept)
            listofPsuedoNumber = range(minOFslope, maxOFslope)
            slopeValue, intercept, r, p, std_erro = linregress(slope_intercept, slopeMassList)


            for everyPsuedo in listofPsuedoNumber:
                medianPoint = everyPsuedo * slopeValue + intercept
                slopeMedian.append(medianPoint)


            takeClosest = lambda num, collection: min(collection, key=lambda x: abs(x - num)) ##### method to find the closest
            firstfilepath = os.path.dirname(file)
            fileName = firstfilepath + "/SlopeFDRfile.txt"
            w = open(fileName, "w")

        #for file in processingFileList:

            with open(file) as PTMfile:
                mod2Xcor = {}
                mainList = ["Scan\tSearchengineRank\tcharge\texpMH\ttheoMH\tExpmz\tXcor\tSeq\tRetentionTime\tProtAcc\tDeltaMod"
                 "\tB-series\tY.series\tIsotpicJump\tdeltaPeptide\tfilename\tCorXcor\tnewExpMH\tLabel\tMedian"
                 "\tCal_Delta_MH\tCal_Delta_M/Z\tCalExp_MZ\t1DawindowCenter\tdecoyCount\ttargetCount\tlocalFDR\t1Da_start\t1Da_end\n"]
                next(PTMfile)
                for newline in PTMfile:
                    if newline!= "\n":
                        splitsnew = newline.split("\t")
                        label = splitsnew[18].strip()
                        delta_modification = float(splitsnew[20].strip())
                        xscore = float(splitsnew[16].strip())

                        if delta_modification > (minOFslope - 0.6) and delta_modification < maxOFslope:

                            closestSlope = takeClosest(delta_modification, slopeMedian)
                            if closestSlope not in mod2Xcor:

                                mod2Xcor[closestSlope] = [[xscore, str(label).strip(), newline.strip(), delta_modification]]
                            else:
                                mod2Xcor[closestSlope].append([xscore, str(label).strip(), newline.strip(), delta_modification])

            for mod in mod2Xcor:
                decoy, target = 0.0, 0.0
                mod2Xcor[mod].sort(key=lambda row: (row[0]), reverse=True)
                MinMaxList = [item[3] for item in mod2Xcor[mod]]
                startingPoint = min(MinMaxList)
                endingPoint = max(MinMaxList)

                for eachline in mod2Xcor[mod]:
                    if eachline[1] == "Decoy":
                        decoy += 1.0
                        eachline.append(decoy * float(alpha))
                        eachline.append(target)
                        try:
                            localFDR = float((decoy * float(alpha)) / target)
                        except ZeroDivisionError:
                            localFDR = "1"
                        eachline.append(str(localFDR).strip())
                        eachline.append(str(startingPoint))
                        eachline.append(str(endingPoint))

                    elif eachline[1] == "Target":
                        target += 1.0
                        eachline.append(decoy * float(alpha))
                        eachline.append(target)
                        try:
                            localFDR = float((decoy * float(alpha)) / target)
                        except ZeroDivisionError:
                            localFDR = "1"
                        eachline.append(str(localFDR))
                        eachline.append(str(startingPoint))
                        eachline.append(str(endingPoint))

            for fdrline in mod2Xcor:
                for everyLine in mod2Xcor[fdrline]:
                    mainList.append("\t".join(everyLine[2].split("\t")) + "\t" + str(fdrline) + "\t" + str(everyLine[4]).strip() + "\t" + str(
                        everyLine[5]).strip() + "\t" + str(everyLine[6]).strip()+ "\t" + str(everyLine[7]).strip() + "\t" + str(everyLine[8]).strip()+ "\n")

            for newline in mainList:
                w.write(str(newline))
            w.close()
