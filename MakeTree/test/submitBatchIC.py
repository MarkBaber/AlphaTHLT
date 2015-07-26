#!/usr/bin/env python

import sys
import os.path
import time
import math
from datetime import date

# ****************************************************************************************************
# *                                               Main                                               *
# ****************************************************************************************************



def main():


    # Configuration parameters
    # ****************************************

    baseDir = "/vols/cms04/mb1512/Batch"
    
    # ****************************************


    # Information to extract
    cfgFile = ""
    jobs    = 0
    suffix  = ""

    argCount = 0

    # Loop through and extract arguments 
    for arg in sys.argv:
        argCount += 1

        if (argCount == 1):         # Skip the python executable call
            continue
        elif (argCount == 2):
            cfgFile = str(os.getcwd()) + "/" + arg
        elif (argCount == 3):
            jobs = int(arg)
        elif (argCount == 4):
            suffix = str(arg)
        else:
            print "Error: More than 3 arguments specified, expected only three: cfgFilename, numberOfJobs, (suffix)"
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

    # Load config file to extract configuration parameters
    import MakeTrees_cfg
    fileCount = len(MakeTrees_cfg.selectedSample.files)
    fileName  = MakeTrees_cfg.selectedSample.name


    # Calculate the required number of files per job
    filesPerJob = int( math.ceil( float(fileCount) / jobs ) )

    # Correct the number of jobs submitted if necessary
    filesToRun  = jobs * filesPerJob
    excessFiles = filesToRun % fileCount
    excessJobs  = int(math.floor( float(excessFiles)/filesPerJob ))

    jobs        -= excessJobs
 
    # Extract CRAB job name
    # ----------------------------------------
    crabJob   = fileName.pythonValue().replace("'","").replace(".root","")


    # Create directory to store jobs
    # ----------------------------------------

    # Directory to return the jobs to
    outputDir = baseDir + "/" + theDate + "_" + crabJob + "_" + suffix

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

    # Make directory
    os.makedirs( outputDir )
    os.makedirs( configOutDir )
                

    # ********************************************************************************
    # *                         Loop over jobs to be created                         *
    # ********************************************************************************


    for jobNum in range(1, jobs + 1):

        # Determine which file range to run over
        indexLow  = (jobNum - 1)*filesPerJob
        indexHigh = (jobNum    )*filesPerJob - 1

        # Create a unique ouput ROOT filename
        outputROOTName = crabJob + "_" + str(jobNum) + ".root"

        # Create a new configuration file for the specific job
        tempCfgText = cfgText

        # Add lowest and highest files to be run over to config file
        header   = ""
        header  += "# <PYTHON JOB STUFF> \n"
        header  += "# Lowest and highest (exclusive) index file to add, indexHigh = -1  =>  No upper bound\n"
        header  += "indexLow  = " + str(indexLow)  + "\n"
        header  += "indexHigh = " + str(indexHigh) + "\n"
        header  += "\n# </PYTHON JOB STUFF> \n\n\n"
        # Add header to config file
        tempCfgText = header + tempCfgText
        
        # Add file range to file loader
        filenamesStr  = "fileNames = cms.untracked.vstring(\n"
        for index in range( indexLow, indexHigh + 1 ):
            if index == fileCount:
                break

            filenamesStr += "\t'" + MakeTrees_cfg.selectedSample.files[index] + "',\n"
            pass
        filenamesStr += ")"

        tempCfgText = tempCfgText.replace("fileNames = selectedSample.files", filenamesStr ) 
        tempCfgText = tempCfgText.replace("fileName = selectedSample.name",  "fileName = cms.string('" + outputROOTName + "')")


        # Create the configuration file and batchshell script for the job
        batchCfgFile = outputDir + "/MakeTree"   + "_" + str(jobNum) + ".py"
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
        tempShText +=     "export CMSSW_PROJECT_SRC=\"SUSY/UCTHLT/CMSSW_7_4_7/src\"\n"
        tempShText +=     "cd /home/hep/mb1512/$CMSSW_PROJECT_SRC\n\n"

        tempShText +=     "# source cms stuff\n"
        tempShText +=     "source /vols/cms/grid/setup.sh\n"
        tempShText +=     "export SCRAM_ARCH=slc6_amd64_gcc481\n"
        tempShText +=     "eval `scramv1 runtime -sh`\n\n"

        tempShText +=     "cd $OUTPUTDIR\n"
        tempShText +=     "cmsRun " + batchCfgFile + "\n\n"
        

        f = open( batchShFile, "w")
        f.write( tempShText )
        f.close()

        # Set the correct permissions for running on the batch
        os.chmod( batchShFile, 0o774 )



        # Job output
        # ----------------------------------------
        print "ROOT output name = ", outputROOTName
        print "File range       = ", indexLow, "-", indexHigh

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
    print "\t\tCrab job    = ", crabJob
    print "\t\tFiles       = ", fileCount
    print "\t\tFilePerJob  = ", filesPerJob
    print "\t\tOutput dir  = ", outputDir


    pass









# Call the main function
if __name__ == "__main__":
            main()


            
