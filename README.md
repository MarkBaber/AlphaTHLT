AlphaTHLT
=========

Working with CMSSW_6_2_5

Installation instructions:


<pre><code>
  mkdir AlphaTHLT
  
  export SCRAM_ARCH=slc6_amd64_gcc481
  cmsrel CMSSW_7_1_8
  cd CMSSW_7_1_8/src
  # Get UCT emulator
  git cms-merge-topic --unsafe cms-l1t-offline:Emulator_1.2
  
  mkdir AlphaTHLT
  cd AlphaTHLT
  cmsenv
  
  git init
  git remote add origin git@github.com:MarkBaber/AlphaTHLT.git
  git fetch origin
  cd $CMSSW_BASE/src
  scram b -j8
</code></pre>
