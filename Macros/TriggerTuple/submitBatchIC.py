#!/usr/bin/env python

import sys
import os.path
import time
import math
from datetime import date

# ****************************************************************************************************
# *                                               Main                                               *
# ****************************************************************************************************

#
#   ./submitBatchIC.py  triggerTuple.C PU40bx25_HCAL3_HPUV_QCD  100 PU40bx25_HCAL3_HPUV_QCD
#


def main():


    # Configuration parameters
    # ****************************************

    baseDir = "/vols/cms04/mb1512/Batch"
    
    # ****************************************


    # Information to extract
    cfgFile = ""
    jobs    = 0
    sample  = ""

    definedSamples = ["test",

                      "PU40bx50_HCAL3_HPUV_QCD",
                      "PU40bx50_HPUV_QCD",
                      "PU40bx25_HCAL3_HPUV_QCD",
                      "PU40bx25_HPUV_QCD",
                                            
                      "SM",
                      "PU20bx25_HCAL3_Signal",
                      "PU20bx25_Signal",


                      "TTBar", 
                      "DYJets",
                      "NuGun",
                      "QCD30to50", 
                      "QCD50to80", 
                      "QCD80to120",
                      "QCD120to170",
                      "QCD170to300",
                      "QCD300to470",
                      "QCD470to600",
                      "QCD600to800",
                      "QCD800to1000",
                      "T2cc_250_210",
                      "T2tt_500_250",
                      "T2tt_300_200",
                      "T1bbbb_2J_mGl_1000_mLSP_900",
                      "T1tttt_2J_mGl_1200_mLSP_800",
                      "T2tt_2J_mStop_425_mLSP_325",
                      "T2tt_2J_mStop_850_mLSP_100",
                      "T2bb_2J_mStop_600_mLSP_580",  
                      "T2tt_2J_mStop_650_mLSP_325",
                      "T2qq_2J_mStop_600_mLSP_550",  
                      "T2tt_2J_mStop_500_mLSP_325",  
                      ]

    argCount = 0

    # Loop through and extract arguments 
    for arg in sys.argv:
        argCount += 1

        if (argCount == 1):         # Skip the python executable call
            continue
        elif (argCount == 2):
            cfgFile = str(os.getcwd()) + "/" + arg
        elif (argCount == 3):
            baseDir += "/" + arg
        elif (argCount == 4):
            jobs = int(arg)
        elif (argCount == 5):
            sample = str(arg)

            if (sample not in definedSamples):
                print "Error: Sample '", sample, "' not recognised, check against the list of approved samples.\n"
                return
        else:
            print "Error: More than 4 arguments specified, expected only four: cfgFilename, batchDirectory, numberOfJobs, sample"
            return
        pass
    
    if ( cfgFile == "" ):
        print "Error: Insufficient arguments specified, expected two: cfgFilename, numberOfJobs"
        return

    print "\n\n"

    # ------------------------------------------------------------------------------------------------------------------------


    # Get current date
    theDate = str(date.today())


    # Open the configuration file
    f = open( cfgFile )
    cfgText = f.read()
    f.close()


    # Create directory to store jobs
    # ----------------------------------------

    # Directory to return the jobs to
    outputDir = baseDir + "/" + theDate + "_" + sample

    # Make sure directory does not exist
    if ( os.path.exists( outputDir ) ):
        # Directory already exists, create a unique filename
        revCount = 0        
        revStr   = ""
        while ( os.path.exists( outputDir + revStr ) ):
            revCount  += 1
            revStr     = "_rev_" + str(revCount)
            pass

        outputDir += revStr
        pass

    # Create output directory for configuration data
    configOutDir = outputDir + "/Config"
    outputSubDir = outputDir + "/output"
    
    # Make directory
    os.makedirs( outputDir )
    os.makedirs( configOutDir )
    os.makedirs( outputSubDir )

    # ********************************************************************************
    # *                         Loop over jobs to be created                         *
    # ********************************************************************************


    for jobNum in range(0, jobs):

#         # Determine which file range to run over
#         indexLow  = (jobNum - 1)*filesPerJob
#         indexHigh = (jobNum    )*filesPerJob - 1

        # Fractional position of jobs through events
        evLowFact     = (jobNum/float(jobs))
        evHighFact    = ( (jobNum + 1)/float(jobs))
        # Offset from fractional position
        evLowOffset   = 1
        evHighOffset  = 1
        if (jobNum == 0):
            evLowOffset = 0
            pass
#         if (jobNum == jobs - 1):
#             evHighOffset = 0
#             pass

       

        # Create a unique ouput ROOT filename
        outputROOTName = sample + "_" + str(jobNum) + ".root"
        # Name of the executable
        exeName = "rootMacro" + "_" + str(jobNum)

        # Create a new configuration file for the specific job
        tempCfgText = cfgText



        # Add lowest and highest files to be run over to config file
        header   = ""
        header  += "// <JOB STUFF> \n"
        header  += "#define RUN_ON_BATCH\n"

        #
        #  limits: low  = nEvents * (JOBN/NJOBS) + 1 ( if JOBN != 0)
        #          high = nEvents * (JOBN + 1/NJOBS)
        
#         header  += "const int NJOBS = " + str(njobs)  + "\n"
#         header  += "const int JOBN  = " + str(index)  + "\n"
        header  += "float eventLowFact  = " + str(evLowFact)  + ";\n"
        header  += "float eventHighFact = " + str(evHighFact) + ";\n"
        header  += "int eventLowOffset  = " + str(evLowOffset)  + ";\n"
        header  += "int eventHighOffset = " + str(evHighOffset) + ";\n"
        header  += "\n// </ JOB STUFF> \n\n\n"



        # Change the output options for the file
        # fileOut = "fOut = new TFile(\"" + outputSubDir + "/" + outputROOTName + "\",\"RECREATE\");"
        # tempCfgText = tempCfgText.replace("fOut = new TFile(", fileOut + "\n//TFile* fOut = new TFile(")
        
        fileOut = "TString fileName = \"" + outputSubDir + "/\" + sample + \"_" + str(jobNum) + ".root\"" 
        tempCfgText = tempCfgText.replace("TString fileName = outdir + sample + \".root\"", fileOut );
        tempCfgText = tempCfgText.replace("selSampleStrs = sampleStrs[\"BATCH\"];","selSampleStrs = sampleStrs[\"" + sample + "\"];")


        # Make the macro execute itself
        footer   = ''
        footer  += "\n\n"
        footer  += "void " + exeName + "(){\n"
        footer  += "    triggerTuple();\n"
        footer  += "}\n\n"
          
        
        # Add header and footer to config file
        tempCfgText = header + tempCfgText + footer

        


        # Create the configuration file and batchshell script for the job
        batchCfgFile = outputDir + "/" + exeName + ".cpp"
        batchShFile  = outputDir + "/SubmitBatchIC" + "_" + str(jobNum) + ".sh"

        f = open( batchCfgFile, "w")
        f.write( tempCfgText )
        f.close()



        tempShText  =     ""
        tempShText +=     "#!/bin/bash\n\n"
        tempShText +=     "# Modified by the python script\n"
        tempShText +=     "export OUTPUTDIR=\"" + outputDir      + "\"\n"
        tempShText +=     "export FILENAME=\""  + outputROOTName + "\"\n\n"

        tempShText +=     "# IC Batch Job Script\n"
        tempShText +=     "export CMSSW_PROJECT_SRC=\"SUSY/UCTHLT/CMSSW_7_4_0_pre9/src\"\n"

        tempShText +=     "cd /home/hep/mb1512/$CMSSW_PROJECT_SRC\n\n"

        tempShText +=     "# source cms stuff\n"
        tempShText +=     "source /vols/cms/grid/setup.sh\n"
        tempShText +=     "export SCRAM_ARCH=slc6_amd64_gcc481\n"
        tempShText +=     "eval `scramv1 runtime -sh`\n\n"


        tempShText +=     "cd $OUTPUTDIR\n"
        tempShText +=     "root -b " + batchCfgFile + "+\n\n"
        


        f = open( batchShFile, "w")
        f.write( tempShText )
        f.close()

        # Set the correct permissions for running on the batch
        os.chmod( batchShFile, 0o774 )



        # Job output
        # ----------------------------------------
        print "ROOT output name       = ", outputROOTName
        print "Event fractional range = ", evLowFact, "-", evHighFact
        print "Event offsets          = ", evLowOffset, "-", evHighOffset

        # # Submit the job
        # os.system( "qsub -o " + configOutDir + " -e " + configOutDir + " -q hepshort.q " + batchShFile )

        # Submit the job
        submitCommand  = "qsub -o " + configOutDir + " -e " + configOutDir # Change stdout/err location
        submitCommand += " -wd " + configOutDir                            # Change return directory
        submitCommand += " -N " + outputROOTName.replace(".root","")       # Rename job
        submitCommand += " -cwd -q hepshort.q " + batchShFile
        os.system( submitCommand )



        pass

    # Output extracted details
    print "\nConfiguration used:\n"
    print "\tConfiguration file = ", cfgFile 
    print "\tNumber of jobs     = ", jobs
    print "\tBase dir           = ", baseDir
    print "\t\tOutput dir  = ", outputDir


    pass









# Call the main function
if __name__ == "__main__":
            main()


            
