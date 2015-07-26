AlphaTHLT
=========

Installation instructions:


```bash
  cmsrel CMSSW_7_4_7
  cd CMSSW_7_4_7/src
  cmsenv
  git cms-addpkg HLTrigger/Configuration
  git cms-merge-topic cms-tsg-storm:hltUpdatesOnTopOf745plusEpsilon_74X
  git cms-checkdeps -A -a
  scram b -j8

  git init
  git clone -b Run2_CMSSW_7_4_7 --single-branch git@github.com:MarkBaber/AlphaTHLT.git
```

Latest frozen menues:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/ConfDB740#Frozen_menus

/frozen/2015/50ns_5e33/v3.4/HLT/V1: same as /dev/CMSSW_7_4_0/50nsGRun/V152
/frozen/2015/25ns14e33/v3.0/HLT/V1: same as /dev/CMSSW_7_4_0/GRun/V103

Checkout menues:
``` bash
hltGetConfiguration /frozen/2015/25ns14e33/v3.0/HLT/V1 --full --offline --mc --unprescale --process HLT2 --globaltag FALL1374_25V4 --l1-emulator 'stage1,gt' --l1Xml L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml  > hlt_frozen_2015_25ns_14e33_v3p0_HLT_V1.py

hltGetConfiguration /frozen/2015/50ns_5e33/v3.4/HLT/V1 --full --offline --mc --unprescale --process HLT2 --globaltag FALL1374_50V0 --l1-emulator 'stage1,gt' --l1Xml L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml --type=50nsGRun > hlt_frozen_2015_50ns_5e33_v3p4_HLT_V1.py
```

Customise HCAL reconstruction method

50ns = HCAL0
25ns = HCAL3

https://twiki.cern.ch/twiki/bin/viewauth/CMS/TriggerStudiesPAGStudyHCALMethod2
