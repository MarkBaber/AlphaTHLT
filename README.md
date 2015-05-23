AlphaTHLT
=========

Working with CMSSW_7_4_2

Installation instructions:


<pre><code>

  cmsrel CMSSW_7_4_2
  cd CMSSW_7_4_2/src
  cmsenv
  git cms-addpkg HLTrigger/Configuration
  scram b -j8

  mkdir AlphaTHLT
  cd AlphaTHLT
  git remote add origin git@github.com:MarkBaber/AlphaTHLT.git

</code></pre>
