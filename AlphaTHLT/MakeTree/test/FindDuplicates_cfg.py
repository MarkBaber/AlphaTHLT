import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing


# Remove duplicate CRAB jobs that share the same job index
def removeDuplicates(inputFiles, searchStr):

    storedJobs     = []
    uniqueFiles    = []
    duplicateFiles = []
    
    for File in inputFiles:

        startIndex = File.find(searchStr) + len(searchStr)
        endIndex   = File.find("_", startIndex)
       
        jobNumber = File[ startIndex : endIndex ]
        
        if jobNumber not in storedJobs:
            uniqueFiles.append( File )
            storedJobs.append( jobNumber )
        else:
            duplicateFiles.append( File )
            #print jobNumber, "is a duplicate"
        pass
    
    if ( len(duplicateFiles) > 0):
        print "Removed", len(duplicateFiles), "duplicates.\n"
        pass


    return uniqueFiles
    
            


from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V1A_PU40bx25_12Nov14_cfi import * # 25ns, HCAL fix


# Samples:
samples = [QCD30to50,  QCD50to80, QCD80to120, QCD120to170, QCD170to300, QCD300to470,  QCD470to600, QCD600to800, QCD800to1000, 
           TTbar, DYJets, NuGun, ]
           #T2cc_250_210, T2tt_500_250, T2tt_300_200 ]
#samples = [T1bbbb_2J_mGl_1000_mLSP_900, T2tt_2J_mStop_850_mLSP_100, T2tt_2J_mStop_425_mLSP_325, T2tt_2J_mStop_500_mLSP_325, T1tttt_2J_mGl_1200_mLSP_800 ]

selectedSample = QCD300to470


searchStr = "hltReRunResults_"


#print selectedSample.files
unique = removeDuplicates( selectedSample.files, "hltReRunResults_" )


print len(unique)
#print unique


# storedJobs     = []
# uniqueFiles    = []
# duplicateFiles = []



# for File in selectedSample.files:
    
    
#     startIndex = File.find(searchStr) + len(searchStr)
#     endIndex   = File.find("_", startIndex)

#     jobNumber = File[ startIndex : endIndex ]

#     if jobNumber not in storedJobs:
        
#         uniqueFiles.append( File )
#         storedJobs.append( jobNumber )
#     else:
#         duplicateFiles.append( File )
#         print jobNumber, "is a duplicate"



# print "Unique files    = ", len(uniqueFiles)
# print "Duplicate files = ", len(duplicateFiles)


    

