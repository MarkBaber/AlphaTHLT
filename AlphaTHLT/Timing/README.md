HLT timing studies
==================

Use the TSG machines vocms003 or vocms004 for timing, to source CMS-SW:
<pre><code>
source /data/release/cmssw/cmsset_default.sh
</code></pre>

Make a release locally for studies:
<pre><code>
/data/user/$USER
</code></pre>

Check no other jobs are running prior to studies and restrict running to a single core:
<pre><code>
taskset -c 7 cmsRun configuration.py
</code></pre>

The timing per event at 100 kHz is equivalent to a single process on a single core on the vocms003 Ivybridge at 162 ms/ev.