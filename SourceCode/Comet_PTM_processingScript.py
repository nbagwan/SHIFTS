__author__ = "nbagwan"

import glob
import os
import pdb
import numpy
import all_stats


separations = "/"
##### this is the sirst step of shifts, which begins with creating the output folder as per the user input #####
##### takes three input; CorrXcor value from user for selecting the best non-modified PSMs, path to raw data from a master file and output folder name #####
##### the method calualtes the machine error from best nonmodfied PSMs (> scoreInput) and calibrated every PSMs and tags the PSms wheter they are Target or Decoys ####

def ProcessingFile(folder, outfoldername, scoreInput):
    allfiles = glob.glob(os.path.join(folder, "*.txt"))

    proteinFastaColumn = 24
    medianDict = {}
    mainList = {}

    mainList1 = ["Scan\tSearchengineRank\tcharge\texpMH\ttheoMH\tExpmz\tXcor\tSeq\tRetentionTime\tProtAcc\tDeltaMod"
                 "\tB-series\tY.series\tIsotpicJump\tdeltaPeptide\tfilename\tCorXcor\tnewExpMH\tLabel\tMedian"
                 "\tCal_Delta_MH\tCal_Delta_M/Z\tCalExp_MZ\n"]

    colListForWriting = [0, 1, 2, 3, 4, 5, 7, 13, 15, 24, 17, 18, 19, 20, 21] # coloumn to be used in output file

    ver2col = [7, 2, 5, 3, 4, 1, 0, 17, 20] # version: [xcor,charge,expM/Z,expMH,theoMH, SE, Scan, DeltaMOD, jump] coloumn to used for difference calculations

    separations = "/"
    outfolder = folder + separations + outfoldername

    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    subname = folder.split("\\")[-1]
    outfilename = outfolder + separations + str(subname) + "_" + "TargetData_Calibration.txt"
    w = open(outfilename, "w")
    logout = outfolder + separations + str(subname) + "_" + "Median_log.txt"
    w1 = open(logout, "w")
    logFile_List = []
    logDic = {}
    file2deltappm = {}

    for file in allfiles:
        with open(file) as f:
            next(f)
            next(f)
            for line in f:
                line = line.strip("\n")
                if line != "\n":
                    splits = line.strip().split("\t")
                    fn = file.split("\\")[-1].strip()
                    scan = str(splits[ver2col[6]].strip("\" "))
                    scan_fn = scan + "_" + fn
                    delta_m = float(splits[ver2col[3]].strip()) - float(splits[ver2col[4]].strip())
                    if -500 < delta_m < 500.0: ##### selecting the delta mass within 500 Da range
                        if str(splits[ver2col[5]]).strip("\" ") == "1":  # and scan_fn not in SE1_dic: selecting only search engine rank one PSMs
                            scanXcor = float(splits[ver2col[0]].strip())
                            scanCharge = float(splits[ver2col[1]].strip())
                            ScanFN = fn
                            scanTheo_MH = float(splits[ver2col[4]].strip())
                            scanExp_MH = float(splits[ver2col[3]].strip())

                            CorrectedXcor = all_stats.correctedXcor(scanXcor, scanExp_MH, scanCharge) #### from the comet Xcor --> an CorrXcor is calculates, takes in account the avarge peptide mass and charge
                            scanDelta_MH = scanExp_MH - scanTheo_MH
                            scanDelta_MZ = scanDelta_MH / scanCharge
                            proteinHITlist = splits[proteinFastaColumn].strip()
                            scan_Sequence = splits[13].strip().upper()



                            Xcorthershold = float(scoreInput)
                            if CorrectedXcor >= Xcorthershold and -0.1 < scanDelta_MZ < 0.1:  ##### new charge condition. search engine rank fixed
                                delta_ppm = (scanExp_MH - scanTheo_MH) / scanTheo_MH * 1000000 ##### calculating the PPM error
                                if "DECOY" not in splits[proteinFastaColumn].strip("\" ") and "INV" not in splits[proteinFastaColumn].strip("\" "): #### checking for target PSMs
                                    if str(ScanFN) not in file2deltappm:
                                        file2deltappm[ScanFN] = [float(delta_ppm)]
                                        logDic[str(ScanFN)] = [delta_ppm]
                                    else:
                                        file2deltappm[ScanFN].append(float(delta_ppm))
                                        logDic[ScanFN].append(delta_ppm)


                            ##### appeding the all the values calcualtes upstairs####
                            ln = []
                            Continue = True
                            if Continue:
                                for each in colListForWriting:
                                    ln += [splits[each].strip("\" ")]
                                ln += [str(ScanFN)]
                                ln += [str(CorrectedXcor)]
                                ln += [str(scanExp_MH)]
                                if "DECOY_" not in splits[proteinFastaColumn].strip("\" ") and "INV_" not in splits[proteinFastaColumn].strip("\" "):
                                    label = "Target"
                                else:
                                    label = "Decoy"

                                ln += [label]
                                if ScanFN not in medianDict:
                                    medianDict[ScanFN] = [ln]
                                else:
                                    medianDict[ScanFN].append(ln)

            ##### in this part a meadian will be calcualtes for every raw file, from the ppm calcualtes abouve ######
            for fileName in medianDict:
                deltam_median = all_stats.medianMethod(file2deltappm[fileName])  ### a simple method to caluclate median of list
                file2deltappm[fileName] = deltam_median


    ###### writing the lof file, contains every file name , number of PSMs used for median caluclation and the value ########
    for j in logDic:
        logFile_List.append(str(j) + "\t" + "no. of points to calculate median are :" + str(
            len(logDic[j])) + "\t" + "Median is: " + str(numpy.median(sorted(logDic[j]))))

        if len(logDic[j]) < 1000:
            print "WARNING: no of point are few for calibration, TRY a lower Xcor"
            print j, "\t", "no. of points to calculate median are :", str(len(logDic[j])), "\t", "Median is: ", str(
                numpy.median(sorted(logDic[j])))
        else:
            print j, "\t", "no. of points to calculate median are :", str(len(logDic[j])), "\t", "Median is: ", str(
                numpy.median(sorted(logDic[j])))


#####this part to use the median calculated abouve and recalibrate every raw file and PSM in them. ######
    for eachFile in medianDict:
        if eachFile not in mainList:
            mainList[eachFile] = []
        for eachLine in medianDict[eachFile]:
            median = file2deltappm[eachFile]
            charge = float(eachLine[2])
            e_mh_plus = float(eachLine[17])
            t_mh_plus = float(eachLine[4])
            deltaModification = float(eachLine[10])

            e_mZ = float(eachLine[17]) / charge
            temp_Delta_mh = e_mh_plus - t_mh_plus
            delta_mz = temp_Delta_mh / charge
            t_mz = t_mh_plus / charge
            # t_mz = e_mZ - (delta_mz)
            calibrated_MZ_ByComet = e_mZ - (median * e_mZ) / 1000000
            calibrated_Delta_MZ = calibrated_MZ_ByComet - t_mz
            calibrated_Delta_MH_ByComet = calibrated_Delta_MZ * charge


            LineVal = eachLine[0:19] + [str(file2deltappm[eachFile]), str(calibrated_Delta_MH_ByComet), str(calibrated_Delta_MZ), str(calibrated_MZ_ByComet)]
            mainList[eachFile].append(LineVal)

        for index in range(len(mainList[eachFile])):
            mainList1.append("\t".join(mainList[eachFile][index]) + "\n")

########## the all the new produced value will be stored in mainList ####3
######### mainList will be used for writing the output of this step. , output of this step will be used as the basis of all other methods and a input tp opther method #####
    for i in mainList1:
        w.write(i)
    w.close()

    for jj in logFile_List:
        w1.write(str(jj) + "\n")
    w1.close()

    return outfilename