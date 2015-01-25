#!/usr/bin/env python 

import sys
import os.path
import time
import math
from datetime import date

# **************************************************************************************************** 
# *                                               Main                                               * 
# **************************************************************************************************** 


nJobs   = 1
samples = ["T2bb_2J_mStop_600_mLSP_580",  
           "T2tt_2J_mStop_425_mLSP_325", 
           "T2tt_2J_mStop_650_mLSP_325",
           "T2qq_2J_mStop_600_mLSP_550",  
           "T2tt_2J_mStop_500_mLSP_325",  
           "T2tt_2J_mStop_850_mLSP_100",
           ]


def main():

    argCount = 0

    # Loop through and extract arguments 
    for arg in sys.argv:
        argCount += 1

        if (argCount == 1):         # Skip the python executable call
            continue
        elif (argCount == 2):
            batchDir =  arg
        else:
            print "Error: More than 1 argument specified, expected only: batchDirectory"
            return
        pass

    
    for sample in samples:
        if sample == "NuGun":
            nJobs = 10
        else:
            nJobs = 1
            pass

        print "\nSubmitting '", sample, "' to the batch in ", nJobs, " jobs."
        os.system( "./submitBatchIC.py makeSUSYHLTAlphaT.C " + batchDir + " " + str(nJobs) + " " + sample )
        print "./submitBatchIC.py makeSUSYHLTAlphaT.C " + batchDir + " " + str(nJobs) + " " + sample
        pass

    return
        


# Call the main function 
if __name__ == "__main__":
            main()
