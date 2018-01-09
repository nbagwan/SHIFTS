import math
import os
import pdb
import numpy

import all_stats

separations = "/"



def MassTagging2Seq_T(histogramFile, apex_thershold, target):
    filepath = os.path.dirname(histogramFile)
    width = []
    mass = 0
    apexList = []
    if target:
        filename = [filepath + "/Target_Peak_for_sequence_Tagging.txt", filepath + "/Target_ApexList_mass.txt"]
    else:
        filename = [filepath + "/decoy_Peak_for_sequence_Tagging.txt", filepath + "/decoy_ApexList_mass.txt"]
    w2 = open(filename[0], "w")
    ApexList = open(filename[1], "w")
    Apex = []
    outputlist1 = [
        ["Bin", "\t", "Frequency", "\t", "Slope1", "\t", "Slope2", "\t", "peak-Width", "\t", "peak-Apex", "\t",
         "Intercept", "\t", "Tag",
         "\n"]]
    with open(histogramFile) as slopeFile:
        next(slopeFile)
        for line4 in slopeFile:
            if line4 != "\n":
                splits4 = line4.split("\t")
                if splits4[4].strip() == "1":  ### width coloumn
                    width.append(line4.strip())
                    if splits4[5].strip() == "1":  ### apex coloumn
                        mass = splits4[0].strip()
                        apexList.append(mass)
                elif splits4[4].strip() == "0":
                    if len(width) != 0:
                        for each in width:
                            # print each.strip(), "\t", mass
                            outputlist1.append([str(each), "\t", str(mass), "\n"])
                    width = []
                    mass = 0

                if splits4[4].strip() == "1" and splits4[5].strip() == "1" and float(splits4[1]) > float(
                        apex_thershold):
                    Apex.append([str(splits4[6].strip()), "\n"])

    for mas in outputlist1:
        w2.writelines(mas)
    w2.close()

    for apex in Apex:
        ApexList.writelines(apex)
    ApexList.close()

    return filename[1]


######################################################### CometPTM #################################################

def taggingMethodForCOMET_PTM(apexfileNameList, calibrationFile, folder, fasta, slopefdrFile, globalfdrfile, sigmafile, fdrTOuse, localFDRs, peakFDRs):

    outputlsit1 = []
    outputlist2 = [
        ["Scan", "\t", "SearchEngineRank", "\t", "Charge", "\t", "Exp Mh", "\t", "Theo mh+", "\t", "Exp MZ", "\t", "Xcor", "\t",
         "Seq", "\t", "RetentionTIme", "\t", "Protein", "\t", "Delta_MOD" ,"\t", "B_series","\t", "Y_series", "\t", "Jumps", "\t",
         "DeltaPeptide", "\t", "FileName", "\t", "CorXcor" "\t", "New_expMH", "\t", "label", "\t",
         "median", "\t", "Calibrated Delta MH", "\t", "Calibrated Delta MZ","\t", "Calibrated EXP MZ",
         "\t", "modifcationMass", "\t", "Sequence-mass", "\t", "fastaDescription", "\t", "SlopeFDR", "\t",
         "startOFslope", "\t", "endOFslope", "\t", "SlopePeak", "\t", "quantification","\t", "ErrorPPM", "\t",
         "classNumber", "\t", "PTMint", "\t", "TagsforPTM-Orphans", "\t", "FinalSeq_Mass","\n"]]

    # only calculated for the target lines in the quantification file
    #MAD = all_stats.madCalculation(calibrationFile, scoreInput=scoreInput, versionofPD=pdversion)
    MAD = all_stats.fullWidthHalfMaximum(sigmafile)
    GlobalFDR = all_stats.globalfdr(globalfdrfile=globalfdrfile)

    print "Finished Median absolute Deviation for ", calibrationFile, "and the mad is:", MAD

    proteinFast = {}
    for line6 in open(fasta):
        if line6 != "\n":
            if line6.startswith(">"):
                if "DECOY_" not in line6 and "INV" not in line6:
                    splits6 = line6.split("|")
                    if splits6[1].strip() not in proteinFast:
                        proteinFast[splits6[1]] = line6.strip()

                        ####### create a list for target and decoy (big database)

    targetApexList = all_stats.createApexList(apexfileNameList[0])
    decoyApexList = all_stats.createApexList(apexfileNameList[0])


    fnn_sc_chDic = {}
    XcorthersholdDic = {}
    XcorthersholdDicWithValue = {}
    with open(slopefdrFile) as fdrFile:
        next(fdrFile)
        for ln in fdrFile:
            if ln != "\n":
                newsplits = ln.split("\t")
                seq = newsplits[7].strip()
                sc = newsplits[0].strip()
                ch = newsplits[2].strip()
                xc = newsplits[16].strip()
                fnn = newsplits[15].strip().split(".")[0]
                fdr = newsplits[26].strip()
                st = newsplits[27].strip()
                ed = newsplits[28].strip()
                slopePeak = newsplits[23].strip()
                peakFDR = newsplits[32].strip()
                st_end = st + "_" + ed
                fnn_sc_ch = fnn + ".raw" + "_" + sc + "_" + ch
                fnn_sc_chDic_seq = fnn + ".raw" + "_" + sc + "_" + ch + "_" + seq
                # if fnn_sc_ch not in fnn_sc_chDic:
                #     fnn_sc_chDic[fnn_sc_ch] = slopePeak + "\t" + fdr + "\t" + st + "\t" + ed + "\t" + st_end + "\t" + xc + "\t" + peakFDR

###################### local fdr check ###########################################################################################################
                if fnn_sc_chDic_seq not in fnn_sc_chDic:
                    fnn_sc_chDic[fnn_sc_chDic_seq] = slopePeak + "\t" + fdr + "\t" + st + "\t" + ed + "\t" + st_end + "\t" + xc + "\t" + peakFDR

                if float(fdr) < float(localFDRs):
                #if float(fdr) < float(fdrTOuse):
                    if st_end not in XcorthersholdDic:
                        XcorthersholdDic[st_end] = [float(xc)]
                    else:
                        XcorthersholdDic[st_end].append(float(xc))

    ########################################################cheching global FDR Xcor thershold for every 1 da range in mass distibution #####################
    for values in XcorthersholdDic:
        PTMxcorthershold = min(XcorthersholdDic[values])
        test1 = XcorthersholdDic[values]
        if float(PTMxcorthershold) > float(GlobalFDR):
            if values not in XcorthersholdDicWithValue:
                XcorthersholdDicWithValue[values] = [PTMxcorthershold, GlobalFDR]
        else:
            #PTMxcorthershold = float(GlobalFDR)
            if values not in XcorthersholdDicWithValue:
                XcorthersholdDicWithValue[values] = [GlobalFDR, GlobalFDR]


    outfileName = folder + "/AllWithSequence-massTag.txt"
    w3 = open(outfileName, "w")

    notassignedFile = folder + "/NotassignedSequences.txt"
    NotassignedW = open(notassignedFile, "w")

    with open(calibrationFile) as PTMfile1:
        next(PTMfile1)
        for line5 in PTMfile1:
            if line5 != "\n":
                splits5 = line5.split("\t")
                scan5 = splits5[0].strip()
                charge = float(splits5[2])
                xcor = float(splits5[16])
                seq2Use = splits5[7].strip()
                fn = splits5[15].split(".")[0]
                uniq5 = scan5 + "_" + fn.strip() + ".raw"
                Delta_bin5 = float(MAD)
                calibrated_delta_MH = float(splits5[20].strip())
                calibrated_delta_MH_trimmed = all_stats.check(calibrated_delta_MH, 3)

                if "|" in splits5[9].strip():
                    proteinlist = splits5[9].split("|")[1]
                else:
                    proteinlist = splits5[9].split("|")[0]

                fn_scan5_charge = fn.strip() + ".raw" + "_" + scan5 + "_" + str(int(charge))
                fn_scan5_charge_seq2Use = fn.strip() + ".raw" + "_" + scan5 + "_" + str(int(charge)) + "_" + seq2Use

                ##3 checking the lists only for the target if small database else both list if large database
                # if smallDatabasePresent=="1" and splits5[15].strip() == "Target": #### old when i was not considering the smallDB


                Continue = True
                if Continue:
                    # checking the label and assigning hte list that has to be check
                    if splits5[18].strip() == "Target":
                        listname = targetApexList
                    else:
                        listname = decoyApexList

                    quantification = "NA"

                    if "DECOY_" not in proteinlist and "INV" not in proteinlist:
                        if proteinlist in proteinFast:
                            fastaDescription = proteinFast[proteinlist]
                    elif "DECOY_" in proteinlist and "INV" in proteinlist:
                        fastaDescription = "NA"


                    takeClosest = lambda num, collection: min(collection, key=lambda x: abs(x - num))
                    if calibrated_delta_MH >= -500 and calibrated_delta_MH <= 500: #and scan5 == "18730":  # and xcor >= 0.179: ## charge considered

                        Error = "NA"
                        classname = "NA"
                        tagOF_PTM_Orphans = "NA"
                        close = takeClosest(calibrated_delta_MH, listname)

                        start = close - (Delta_bin5 * charge)
                        end = close + (Delta_bin5 * charge)
                        if calibrated_delta_MH >= start and calibrated_delta_MH <= end:
                            #if fn_scan5_charge in fnn_sc_chDic:
                            if fn_scan5_charge_seq2Use in fnn_sc_chDic:
                                start_end = fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[4]
                                fdrvalueforPSM = fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[1]
                                fdrvalueforPEAK = fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[6]

                                if start_end.strip() in XcorthersholdDicWithValue and float(fdrvalueforPSM) < float(localFDRs): #### chaeck for local FDR at 1 da
                                    if xcor >= float(XcorthersholdDicWithValue[start_end][0]): #####check for global FDR
                                        thrs = XcorthersholdDicWithValue[start_end][0]
                                        thrs1 = XcorthersholdDicWithValue[start_end][1]
                                        deltaPeptide = splits5[14].strip()
                                        if deltaPeptide.isalpha():
                                            finalSeqMassTag = deltaPeptide
                                        else:
                                            for index in range(0,len(deltaPeptide)):
                                                if deltaPeptide[index] == "[":
                                                    startofMOD= index
                                                if deltaPeptide[index] == "]":
                                                    endofMOD = index
                                            tempPepe= deltaPeptide[startofMOD+1:endofMOD]
                                            newclose = all_stats.check(close,6)
                                            finalSeqMassTag = deltaPeptide.replace(tempPepe, str(newclose))

                                        Seq_mass = splits5[7].strip() + "_" + str(close)

                                        outputlist2.append(
                                            ["\t".join(line5.strip().split("\t")[:23]), "\t",
                                             "\t".join(line5.strip().split("\t")[30:-1]), str(close), "\t", Seq_mass,
                                             "\t",
                                             fastaDescription,
                                             "\t", fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[1], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[2], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[3], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[0], "\t", quantification.strip(),
                                             "\t", Error, "\t", classname, "\t", str(close), "\t", tagOF_PTM_Orphans,"\t",finalSeqMassTag,
                                             "\n"])

                                elif start_end.strip() in XcorthersholdDicWithValue and float(fdrvalueforPSM) > float(fdrTOuse) and float(fdrvalueforPEAK) < float(peakFDRs): ### check for peak FDR

                                    if xcor >= float(XcorthersholdDicWithValue[start_end][1]): #####check for global FDR

                                        thrs = XcorthersholdDicWithValue[start_end][0]
                                        thrs1 = XcorthersholdDicWithValue[start_end][1]
                                        deltaPeptide = splits5[14].strip()
                                        if deltaPeptide.isalpha():
                                            finalSeqMassTag = deltaPeptide
                                        else:
                                            for index in range(0,len(deltaPeptide)):
                                                if deltaPeptide[index] == "[":
                                                    startofMOD= index
                                                if deltaPeptide[index] == "]":
                                                    endofMOD = index
                                            tempPepe= deltaPeptide[startofMOD+1:endofMOD]
                                            newclose = all_stats.check(close, 6)
                                            finalSeqMassTag = deltaPeptide.replace(tempPepe, str(newclose))
                                        # temp = deltaPeptide.find("[")
                                        #
                                        # if temp!= -1:
                                        #     temp2 = deltaPeptide.find("]", start=temp)
                                        Seq_mass = splits5[7].strip() + "_" + str(close)
                                        # print line5.strip()
                                        # print line5.strip().split("\t")[:23]
                                        # print line5.strip().split("\t")[30:-1]


                                        outputlist2.append(
                                            ["\t".join(line5.strip().split("\t")[:23]), "\t",
                                             "\t".join(line5.strip().split("\t")[30:-1]), str(close), "\t", Seq_mass,
                                             "\t",
                                             fastaDescription,
                                             "\t", fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[1], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[2], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[3], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[0], "\t", quantification.strip(),
                                             "\t", Error, "\t", classname, "\t", str(close), "\t", tagOF_PTM_Orphans,"\t",finalSeqMassTag,
                                             "\n"])
                        else:
                            if fn_scan5_charge_seq2Use in fnn_sc_chDic:
                                start_end = fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[4]
                                fdrvalueforPSM = fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[1]
                                fdrvalueforPEAK = fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[6]
                                if start_end.strip() in XcorthersholdDicWithValue and float(fdrvalueforPSM) < float(localFDRs): #or float(fdrvalueforPEAK) < 0.01:
                                    if xcor >= float(XcorthersholdDicWithValue[start_end][0]):
                                        thrs = XcorthersholdDicWithValue[start_end][0]
                                        thrs1 = XcorthersholdDicWithValue[start_end][1]
                                        Seq_mass = splits5[7].strip() + "_" + str(close)
                                        outputlsit1.append(
                                            ["\t".join(line5.strip().split("\t")[:23]), "\t",
                                             "\t".join(line5.strip().split("\t")[30:-1]), str(close), "\t", Seq_mass,
                                             "\t",
                                             fastaDescription,
                                             "\t", fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[1], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[2], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[3], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[0], "\t", quantification.strip(),
                                             "\n"])

                                elif start_end.strip() in XcorthersholdDicWithValue and float(fdrvalueforPSM) > float(fdrTOuse) and float(fdrvalueforPEAK) < float(peakFDRs):
                                    if xcor >= float(XcorthersholdDicWithValue[start_end][1]):
                                        thrs = XcorthersholdDicWithValue[start_end][0]
                                        thrs1 = XcorthersholdDicWithValue[start_end][1]
                                        Seq_mass = splits5[7].strip() + "_" + str(close)
                                        outputlsit1.append(
                                            ["\t".join(line5.strip().split("\t")[:23]), "\t",
                                             "\t".join(line5.strip().split("\t")[30:-1]), str(close), "\t", Seq_mass,
                                             "\t",
                                             fastaDescription,
                                             "\t", fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[1], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[2], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[3], "\t",
                                             fnn_sc_chDic[fn_scan5_charge_seq2Use].split("\t")[0], "\t", quantification.strip(),
                                             "\n"])

    for tag in outputlist2:
        w3.writelines(tag)
    w3.close()

    # return outfileName

    for notassigned in outputlsit1:
        NotassignedW.writelines(notassigned)
    NotassignedW.close()