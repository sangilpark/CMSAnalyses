if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./test-run-9x-mc.sh EVENTS"
    echo "  "
    exit -1
fi

export EVENTS=$1


#
# GlobalTag choice
#
# cat /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_5/src/Configuration/AlCa/python/autoCond.py | grep run2_mc
#    'run2_mc_50ns'      :   '80X_mcRun2_startup_v12',
#    'run2_mc'           :   '80X_mcRun2_asymptotic_v12',
#


export MYFILE=file:/u/user/sangilpark/WorkDir/CP_Violation/TestSample/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1_MINIAODSIM.root
### Type              Module                  Label   Process
### LHEEventProduct   "externalLHEProducer"   ""      "LHE"


# python ssbanalyzer_cfg.py print                   \
cmsRun ssbanalyzer_cfg.py print                   \
    label=DYJetsToLL			 \
    outputFile=SSBTree.root    		 \
    globalTag=94X_mc2017_realistic_v14	 \
    maxEvents=${EVENTS}        		 \
    inputFiles=${MYFILE}

