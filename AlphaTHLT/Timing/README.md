HLT timing studies
==================

Use the TSG machines vocms003 or vocms004 for timing, to source CMS-SW:
<code>
source /data/release/cmssw/cmsset_default.sh
</code>

Make a release locally for studies:
<code>
/data/user/$USER
</code>

Check no other jobs are running prior to studies and restrict running to a single core:
<code>
taskset -c 7 cmsRun configuration.py
</code>

The timing per event at 100 kHz is equivalent to a single process on a single core on the vocms003 Ivybridge at 162 ms/ev.

Further reading: 

    https://twiki.cern.ch/twiki/bin/viewauth/CMS/TriggerStudiesTiming
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/FastTimerService