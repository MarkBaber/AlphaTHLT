AlphaTHLT
=========

Working with CMSSW_7_4_2

Installation instructions:


```bash
  cmsrel CMSSW_7_4_5
  cd CMSSW_7_4_5/src
  cmsenv
  git cms-addpkg HLTrigger/Configuration
  git cms-merge-topic cms-tsg-storm:hltUpdatesOnTopOf744_74X
  git cms-checkdeps -A -a
  scram b -j8

  git init
  git clone git@github.com:MarkBaber/AlphaTHLT.git
```

Latest frozen menues:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/ConfDB740#Frozen_menus

/frozen/2015/50ns_5e33/v2.1/HLT/V5  same as /dev/CMSSW_7_4_0/50nsGRun/V120 (but customized for the online)
/frozen/2015/25ns14e33/v2.1/HLT/V1  same as /dev/CMSSW_7_4_0/GRun/V76


Checkout menues:
``` bash
hltGetConfiguration /frozen/2015/25ns14e33/v2.1/HLT/V1 --full --offline --mc --unprescale --process HLT2 --globaltag FALL1374_25V4 --l1-emulator 'stage1,gt' --l1Xml L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml  > hlt_frozen_2015_25ns_14e33_v2p0_HLT_V1.py

hltGetConfiguration /frozen/2015/50ns_5e33/v2.1/HLT/V5 --full --offline --mc --unprescale --process HLT2 --globaltag FALL1374_50V0 --l1-emulator 'stage1,gt' --l1Xml L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml --type=50nsGRun > hlt_frozen_2015_50ns_5e33_v2p0_HLT_V1.py
```

Customise HCAL reconstruction method

50ns = HCAL0
25ns = HCAL3

https://twiki.cern.ch/twiki/bin/viewauth/CMS/TriggerStudiesPAGStudyHCALMethod2
