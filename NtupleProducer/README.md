Basic Instructions

```
cmsrel CMSSW_9_2_0
cd CMSSW_9_2_0/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v1.14.1
git clone git@github.com:jngadiub/FastPUPPI.git -b pf-ml-l1
scram b -j8
```

To produce the L1PF inputs:
```
cd FastPUPPI/NtupleProducer/python/
cmsRun FastPUPPI/NtupleProducer/python/runInputs.py
```

Example how to submit jobs to the lxbatch:
```
cd FastPUPPI/NtupleProducer/python
./scripts/cmsSplit.pl --dataset /SinglePion0_FlatPt-8to100/PhaseIISpring17D-PU140_90X_upgrade2023_realistic_v9-v1/GEN-SIM-DIGI-RAW --files-per-job 1 --label SinglePion0_PU140 runInputs.py --lsf 8nh --eosoutdir /eos/cms/store/cmst3/group/dehep/L1PFInputs/SinglePion0_PU140
runInputs_SinglePion0_PU140_bsub.sh
```
For more options for job splitting:
```
./scripts/cmsSplit.pl --help
```
The job outputs are copied to the eos directory --eosoutdir

When jobs are finished you can clean up your local dir
```
./runInputs_SinglePion0_PU140_cleanup.sh
```

The trigger MC can be found on DAS `dataset=/*/*PhaseIISpring17D*/*`

Other resources: <br>
[Correlator code](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_9_2_0_and_l1t_phase2_v1_10) <br>
[Track trigger code](https://twiki.cern.ch/twiki/bin/view/CMS/L1Tracklet90X#Recipe_for_CMSSW_9_2_0) <br>
