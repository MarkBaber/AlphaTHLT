AlphaTHLT
=========

Working with CMSSW_7_4_2

Installation instructions:


<pre><code>

  cmsrel CMSSW_7_4_2
  cd CMSSW_7_4_2/src
  cmsenv
  git cms-addpkg HLTrigger/Configuration
  git cms-merge-topic cms-tsg-storm:hltUpdatesOnTopOf742
  scram b -j8

  git init
  git clone git@github.com:MarkBaber/AlphaTHLT.git

</code></pre>


/frozen/2015/50ns_5e33/v2.0/HLT/V1  same as /dev/CMSSW_7_4_0/50nsGRun/V92
/frozen/2015/25ns_14e33/v2.0/HLT/V1 same as /dev/CMSSW_7_4_0/GRun/V60

Checkout menues:
https://twiki.cern.ch/twiki/bin/viewauth/CMS/ConfDB740#Frozen_menus

<pre><code>


# hltGetConfiguration /frozen/2015/25ns_14e33/v2.0/HLT/V1 --full --offline --mc --unprescale --process HLT2 --globaltag FALL1374_25V4 --l1-emulator 'stage1,gt' --l1Xml L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml  > hlt_frozen_2015_25ns_14e33_v2p0_HLT_V1.py

hltGetConfiguration /frozen/2015/25ns14e33/v2.0/HLT/V2 --full --offline --mc --unprescale --process HLT2 --globaltag FALL1374_25V4 --l1-emulator 'stage1,gt' --l1Xml L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml  > hlt_frozen_2015_25ns_14e33_v2p0_HLT_V1.py


hltGetConfiguration /frozen/2015/50ns_5e33/v2.0/HLT/V3 --full --offline --mc --unprescale --process HLT2 --globaltag FALL1374_50V0 --l1-emulator 'stage1,gt' --l1Xml L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml --type=50nsGRun > hlt_frozen_2015_50ns_5e33_v2p0_HLT_V1.py




</code></pre>


50ns = HCAL0
25ns = HCAL3

https://twiki.cern.ch/twiki/bin/viewauth/CMS/TriggerStudiesPAGStudyHCALMethod2
