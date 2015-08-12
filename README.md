AlphaTHLT
=========

Installation instructions:


```bash
  cmsrel CMSSW_7_4_8_patch1
  cd CMSSW_7_4_8_patch1/src
  cmsenv
  git cms-addpkg HLTrigger/Configuration
  scram build
  cd HLTrigger/Configuration/test
  ./cmsDriver.csh   

  git init
  git clone -b Run2_CMSSW_7_4_8_patch1 --single-branch git@github.com:MarkBaber/AlphaTHLT.git
```

Data: Taken from HLT_50ns_5e33_v3_cff.py


Latest frozen menues:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/ConfDB740#Frozen_menus

