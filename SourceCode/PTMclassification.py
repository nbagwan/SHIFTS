__author__ = "nbagwan"
import pdb
import all_stats
import numpy


def ptmClassifier(listofpaths, sigmaFIle):
    mainList = []
    finalList = []
    for everypath in listofpaths:
        filename = everypath + "/" + "NotassignedSequences.txt"
        with open(filename) as unassignedFile:
            for line in unassignedFile:
                if line != "\n":
                    splits = line.split("\t")
                    calibratedDelta_MH = splits[20].strip()
                    mainList.append([everypath + "\t" + line.strip() + "\t", float(calibratedDelta_MH.strip())])

    mainList.sort(key=lambda x: x[1])
    for index in range(len(mainList)):
        m1 = mainList[index][1]
        if index == 0:
            m0 = 0
        else:
            m0 = mainList[index - 1][1]
        try:
            ppmError = abs((m0 - m1) / m0) * 1000000
        except ZeroDivisionError:
            ppmError = 0.0

        finalList.append([mainList[index][0], str(ppmError)])

    SigmaValue = all_stats.sigmaFinder(listofpaths[0] + "/" + sigmaFIle)
    sigmaRange = float(SigmaValue)

    classList = []
    for index1 in range(len(finalList)):
        if index1 == 0:
            classificationTerm = 1
            classList.append([finalList[index1][0] + finalList[index1][1], str(classificationTerm)])
        else:
            if float(finalList[index1][1]) <= sigmaRange:
                classificationTerm = classificationTerm
                classList.append([finalList[index1][0] + finalList[index1][1], str(classificationTerm)])
            else:
                classificationTerm = classificationTerm + 1
                classList.append([finalList[index1][0] + finalList[index1][1], str(classificationTerm)])

    classDic = {}
    for classes in classList:
        if classes[1].strip() not in classDic:
            classDic[classes[1].strip()] = [float(classes[0].split("\t")[21])]
        else:
            classDic[classes[1].strip()].append(float(classes[0].split("\t")[21]))


    for dikey in classDic:
        medianValuelist = classDic[dikey]
        if len(medianValuelist) < 2:
            massOfOrphan = medianValuelist[0]
            #intformula = int(massOfOrphan / float(binvalue)) * float(binvalue) + float(binvalue) / 2

            orphanPTM = all_stats.check(massOfOrphan, 6)
            classDic[dikey].append(orphanPTM)
        else:
            massOfOrphan = numpy.median(medianValuelist)
            #intformula = int(massOfOrphan / float(binvalue)) * float(binvalue) + float(binvalue) / 2
            orphanPTM = all_stats.check(massOfOrphan, 6)
            classDic[dikey].append(orphanPTM)

    for appendingPTMinList in classList:
        if appendingPTMinList[1] in classDic:
            ptmINT = classDic[appendingPTMinList[1]][1:][0]
            appendingPTMinList.append(all_stats.check(float(ptmINT),6))
            appendingPTMinList.append("Orphan")

    DicTomergeFiles = {}
    for ii in classList:
        PathTO_ptmFile = ii[0].split("\t")[0]
        if PathTO_ptmFile not in DicTomergeFiles:
            DicTomergeFiles[PathTO_ptmFile] = [ii]
        else:
            DicTomergeFiles[PathTO_ptmFile].append(ii)

    for pathInDic in DicTomergeFiles:
        if pathInDic in listofpaths:
            tag = pathInDic + "/" + "AllWithSequence-massTag.txt"
            tagFile = open(tag, "a")
            for linestring in DicTomergeFiles[pathInDic]:
                deltaPeptide = linestring[0].split("\t")[15]
                if deltaPeptide.isalpha():
                    finalSeqMassTag = deltaPeptide
                else:
                    for index in range(0, len(deltaPeptide)):
                        if deltaPeptide[index] == "[":
                            startofMOD = index
                        if deltaPeptide[index] == "]":
                            endofMOD = index
                    tempPepe = deltaPeptide[startofMOD + 1:endofMOD]
                    newclose = all_stats.check(float(linestring[2]), 6)
                    finalSeqMassTag = deltaPeptide.replace(tempPepe, str(newclose))


                lineToadd = "\t".join(linestring[0].split("\t")[1:]) + "\t" + str(linestring[1]) + "\t" + str(linestring[2]) + "\t" + str(linestring[3]) + "\t" + finalSeqMassTag + "\n"
                tagFile.writelines(lineToadd)

            tagFile.close()

#ptmClassifier(listofpaths=['E:\\Users\\nbagwan\\Desktop\\OpenSearchPaper\\Improved_isPTM\\tmt-1/Mad_testHW_global',
   #                        'E:\\Users\\nbagwan\\Desktop\\OpenSearchPaper\\Improved_isPTM\\tmt-2/Mad_testHW_global'],sigmaFIle="sigmaCalculations.txt", binvalue= 0.001)