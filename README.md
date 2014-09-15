AlphaTHLT
=========

Working with CMSSW_6_2_5

Installation instructions:


<pre><code>
  mkdir AlphaTHLT
  cmsrel CMSSW_X_Y_Z
  cd CMSSW_X_Y_Z/src
  mkdir AlphaTHLT
  cmsenv
  
  git init
  git remote add origin git@github.com:MarkBaber/AlphaTHLT.git
  git pull origin master
  scram b -j8
</code></pre>
