AlphaTHLT
=========

  mkdir AlphaTHLT
  cd AlphaTHLT
  cmsenv
  
  git init
  git remote add origin git@github.com:MarkBaber/AlphaTHLT.git
  git fetch origin
  cd $CMSSW_BASE/src
  scram b -j8
