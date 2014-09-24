AlphaTHLT
=========

Working with CMSSW_7_1_0_pre8

Installation instructions:


<pre><code>
  mkdir AlphaTHLT
  cmsrel CMSSW_7_1_0_pre8
  cd CMSSW_7_1_0_pre8/src
  cmsenv
  git cms-merge-topic --unsafe cms-l1t-offline:Emulator_1.2

  mkdir AlphaTHLT
  cd AlphaTHLT
  git init
  git remote add origin git@github.com:MarkBaber/AlphaTHLT.git
  git fetch origin
  git checkout Run2_CMSSW_7_1_0_pre8
  cd $CMSSW_BASE/src
  scram b -j8
</code></pre>
