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

1. HLT emulation
=========
The first stage for trigger measurements is the remulation of the HLT using from CMSSW RAW datasets. This utilises configurations in `AlphaTHLT/ReRunHLT/test/` which are generated using the `hltGetConfiguration` command described in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT.

Running jobs
------------
You can submit grid jobs using the CRAB configuration files in: `AlphaTHLT/ReRunHLT/crab/`. This will submit GRID jobs to process the trigger emulation, returning the processed files to Imperial, edit as required to return to your area on the Imperial Tier-2. The name of the basedirectory of the CRAB jobs is used as the key for the next step of processing. You can submit your chosen configuration uisng the following example command:
  ```
  crab submit crab3_Data_Run2015D.py
  ```

2. Making trigger Ntuples
=========
Trigger flat trees are made from processing of CMSSW Ntuples produced in the previous step.

Creating sample files for processed files
--------------------
At Imperial setup the grid environment and generate a GRID proxy in order to access the Tier-2:
  ```
  source /vols/cms/grid/setup.sh
  voms-proxy-init --voms cms --valid 168:00
  ```
Then execute `getCRAB3Jobs.py` using the name of the CRAB basedirectory:
  ```
  /home/hep/mb1512/.scripts/Jobs/getCRAB3Jobs.py <NAME_OF_CRABDIR>
  ```
This will automatically search through directories on the Tier-2 for all files associated with the CRAB job, producing filtered CMSSW `cfi` files with links to all the processed files stored in a PSet. After generation perform a `scram build`:
```
  cd ../../..
  scram b -j8
```
Creating trigger Ntuples
--------------------
Trigger ntuples are produced by the configuration file `MakeTrees_cfg.py` in the directory `AlphaTHLT/MakeTree/test`. To run interactively use the command:
```
cmsRun MakeTrees_cfg.py
```

You can browse files on the Tier-2 with the command:
  ```
  lcg-ls $DCACHE_SRM_ROOT/store/user/
  ```
  Remove files with the command:
  ```
  lcg-del -l FILEPATH
   ```
   And remove empty directories with:
    ```
  lcg-del -a FILEPATH
   ```
