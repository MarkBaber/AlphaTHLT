AlphaTHLT
=========

Installation instructions:


```bash
cmsrel CMSSW_8_0_3_patch1
cd CMSSW_8_0_3_patch1/src
cmsenv
rehash

git cms-init
git cms-addpkg HLTrigger/Configuration
git cms-addpkg L1Trigger/L1TCalorimeter L1Trigger/L1TMuon L1Trigger/L1TGlobal
git cms-merge-topic -u 13704
git cms-merge-topic -u 13767
git cms-merge-topic -u cms-tsg-storm:HLTinApril80X
git cms-merge-topic -u 13809
git cms-merge-topic -u cms-tsg-storm:L1T_externals_for_803

git cms-checkdeps -A -a
scram b -j 4
cmsenv
rehash

  git init
  git clone -b Run2_CMSSW_8_0_3_patch1 --single-branch git@github.com:MarkBaber/AlphaTHLT.git
```
