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
samples = ["TTBar",
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
           "QCD800to1000"
           ]


def main():
    
    for sample in samples:
        if sample == "NuGun":
            nJobs = 10
        else:
            nJobs = 1
            pass

        print "\nSubmitting '", sample, "' to the batch in ", nJobs, " jobs."
        os.system( "./submitBatchIC.py makeSUSYHLTAlphaT.C " + str(nJobs) + " " + sample )
        pass

    return
        


# Call the main function 
if __name__ == "__main__":
            main()
