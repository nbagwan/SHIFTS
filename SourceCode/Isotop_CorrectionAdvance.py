__author__ = "nbagwan"

import math
import pdb
import all_stats
import re

# ReadFile = "TargetData_withSequence-massTag.txt"
# outputFile = "/IsotopCorrection_TargetData_withSequence-massTag.txt"


# ReadFile="inputfile_mass.txt"
# outputFile="outputfile_mass.txt"

def separateTargetAndDecoy(fileName, folder):
    Tdic = {}
    Ddic = {}
    filetouse = folder + "/" + fileName
    with open(filetouse) as filename:
        next(filename)
        for line in filename:
            splits = line.split("\t")
            label = splits[18].strip()

            if label == "Target":
                Tdic[line.strip()] = 1
            else:
                Ddic[line.strip()] = 1

    return Tdic, Ddic


def createDic(Dic, key, val, header):
    neg_dic = {}
    pos_dic = {}
    for line in Dic:
        splits = line.split("\t")
        name = splits[key].strip()
        value = float(splits[val].strip())

        if name not in neg_dic:
            neg_dic[name] = [[float(splits[val].strip()), line.strip()]]
    #         pdb.set_trace()
        else:
            neg_dic[name].append([float(splits[val].strip()), line.strip()])
        #if value < 0:

        #     if name not in neg_dic:
        #         neg_dic[name] = [[float(splits[val].strip()), line.strip()]]
        #         pdb.set_trace()
        #     else:
        #         neg_dic[name].append([float(splits[val].strip()), line.strip()])
        # else:
        #     if name not in pos_dic:
        #         pos_dic[name] = [[float(splits[val].strip()), line.strip()]]
        #     else:
        #         pos_dic[name].append([float(splits[val].strip()), line.strip()])

    return neg_dic


def checkDifference(seq2massDic, MainList, FWHMvalue):
    c13Mass = 1.003354
    c13end = c13Mass + float(FWHMvalue)
    c13start = c13Mass - float(FWHMvalue)

    sulphurMass = 0.9934
    sulphurEnd = sulphurMass + float(FWHMvalue)
    sulphurStart = sulphurMass - float(FWHMvalue)

    for eachSeq in seq2massDic:
        AllMassList = seq2massDic[eachSeq]
        AllMassList.sort(key=lambda x: x[0])  ## sorting the list with mass for each sequence
        massList = [item[0] for item in AllMassList]
        LineList = [item[1] for item in AllMassList]
        tagList = ["no"] * len(AllMassList)
        updateAllMassList = map(lambda x, y, z: [x, y, z], massList, LineList, tagList)
        if checkingDiffMass(massList) == len(massList):
            MainList.append(updateAllMassList)
        else:
            newList = [massList[0]]
            tagList = ["no"]
            counter = 0
            sulphur = False
            for index in range(len(AllMassList) - 1):
                if sulphur:
                    mass1 = massList[index - 1]
                else:
                    mass1 = massList[index]
                mass2 = massList[index + 1]
                difference = round(mass2 - mass1, 4)
                # difference = mass2 - mass1
                # print difference, mass1, mass2, index, len(newList),newList[index]
                # pdb.set_trace()
                if counter < 2:
                    if c13end >= float(difference) >= c13start:
                        # newList[index]=mass2
                        counter += 1
                        mass2 = newList[index]
                        if checkingDiffMass(massList[len(newList):len(massList)]) == len(massList) - len(newList):
                            newList += [mass2] * (len(massList) - len(newList))
                            tagList += ["yes"] * (len(massList) - len(tagList))
                            # pdb.set_trace()
                            # MainList.append(map(lambda x,y:[x,y],newList,LineList))
                            break
                        else:
                            newList.append(mass2)
                            tagList.append("yes")
                    elif difference == 0 and newList[index] != mass1:
                        newList.append(newList[index])
                        tagList.append("yes")
                    elif counter == 1 and sulphurEnd >= float(difference) >= sulphurStart and newList[index] == newList[
                                index - 1] and \
                                    newList[index] != mass1:
                        tagList.append("sulfur")
                        sulphur = True
                        newList.append(newList[index])
                    else:
                        newList.append(mass2)
                        tagList.append("no")
                else:
                    newList.append(mass2)
                    tagList.append("no")
            # if len(tagList)!=len(newList):
            #     print newList, tagList
            MainList.append(map(lambda x, y, z: [x, y, z], newList, LineList, tagList))
    # print MainList
    return MainList


def checkingDiffMass(list):
    count = 0
    old = 10000
    for each in list:
        if old == 10000:
            old = each
        else:
            new = each
            if new == old:
                count += 1
            old = new
    return count + 1


def isotop(targetFile, folder, madFIle):
    header = ["Scan", "SearchEngineRank", "Charge", "Exp Mh", "Theo mh+", "Exp MZ", "Xcor",
                      "Seq", "RetentionTIme", "Protein", "DeltaMod","B_series","Y_series","Jumps","DeltaPeptide", "FileName", "CorXcor", "new_expMH","label",
                      "median", "Calibrated Delta MH","Calibrated Delta MZ","Calibrated EXP MZ",
                      "modifcationMass", "Sequence-mass", "fastaDescription", "SlopeFDR",
                      "startOFslope", "endOFslope", "SlopePeak", "quantification",
                      "ErrorPPM",
                      "classNumber", "PTMint", "TagsforPTM-Orphans", "FinalSeq_Mass",
                      "Corr_Seq-mass", "Corr_mass", "Monoisotop_T/F"]

    writefile = folder + "/IsotopCorrection_TargetData_withSequence-massTag.txt"
    wfile = open(writefile, "w")
    wfile.write("\t".join(header) + "\n")

    FWHM = float(all_stats.fullWidthHalfMaximum(MADfile=madFIle))

    Tdic, Ddic = separateTargetAndDecoy(targetFile, folder=folder)


    Seq2MassNegDic = createDic(Tdic, key=7, val= 33, header =True)

    MainList = []
    #MainList = checkDifference(Seq2MassPosDic, MainList, FWHM)
    MainList = checkDifference(Seq2MassNegDic, MainList, FWHM)
    for lists in MainList:
        for list2 in lists:
            deltaPeptide = list2[1].split("\t")[35]
            cometOfficialSeq_teg = list2[1].split("\t")[7]
            newclose = all_stats.check(list2[0], 6)
            if "[" in deltaPeptide:
                for index in range(0, len(deltaPeptide)):
                    if deltaPeptide[index] == "[":
                        startofMOD = index
                    if deltaPeptide[index] == "]":
                        endofMOD = index
                tempPepe = deltaPeptide[startofMOD + 1:endofMOD]
                #newclose = all_stats.check(list2[0], 2)
                finalSeqMassTag = deltaPeptide.replace(tempPepe, str(newclose))

            else:
                finalSeqMassTag = cometOfficialSeq_teg + "_" + str(newclose)
            #wfile.write(list2[1] + "\t" + list2[1].split("\t")[7] + "_" + str(list2[0]) + "\t" + str(list2[0]) + "\t" +
            #list2[2] + "\n")

            wfile.write(
                list2[1] + "\t" + finalSeqMassTag + "\t" + str(list2[0]) + "\t" +
                list2[2] + "\n")



     # for decoy
    Seq2MassNegDic = createDic(Ddic, key=7, val=33,
                                   header=True)
    MainList = []
    #MainList = checkDifference(Seq2MassPosDic, MainList, FWHM)
    MainList = checkDifference(Seq2MassNegDic, MainList, FWHM)

    for lists in MainList:
        for list2 in lists:
                #wfile.write(list2[1] + "\t" + list2[1].split("\t")[7] + "_" + str(list2[0]) + "\t" + str(list2[0]) + "\t" +
                #list2[2] + "\n")
            deltaPeptide = list2[1].split("\t")[35]
            cometOfficialSeq_teg = list2[1].split("\t")[7]
            newclose = all_stats.check(list2[0], 6)
            if "[" in deltaPeptide:
                for index in range(0, len(deltaPeptide)):
                    if deltaPeptide[index] == "[":
                        startofMOD = index
                    if deltaPeptide[index] == "]":
                        endofMOD = index
                tempPepe = deltaPeptide[startofMOD + 1:endofMOD]
                finalSeqMassTag = deltaPeptide.replace(tempPepe, str(newclose))

            else:
                finalSeqMassTag = cometOfficialSeq_teg + "_" + str(newclose)
            # wfile.write(list2[1] + "\t" + list2[1].split("\t")[7] + "_" + str(list2[0]) + "\t" + str(list2[0]) + "\t" +
            # list2[2] + "\n")

            wfile.write(
                list2[1] + "\t" + finalSeqMassTag + "\t" + str(list2[0]) + "\t" +
                list2[2] + "\n")


    wfile.close()
    return writefile


# isotop(targetFile="AllWithSequence-massTag.txt",
#        folder="E:/Users/nbagwan/Desktop/OpenSearchPaper/Improved_isPTM/tmt-1/Mad_testHW_global_all",
#        experimentname="TMT10", smallDbP="1",
#        pdversion="2.1")
