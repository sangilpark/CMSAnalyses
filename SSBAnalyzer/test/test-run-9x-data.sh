if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./test-run-9x-data.sh EVENTS"
    echo "  "
    exit -1
fi

export EVENTS=$1


#
# GlobalTag choice
#
# cat /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_17/src/Configuration/AlCa/python/autoCond.py | grep run2_data
#    'run2_data'         :   '80X_dataRun2_v17',
# 


export MYFILE=file:////u/user/sangilpark/WorkDir/CP_Violation/TestSample/DoubleMuon_Run2017B-31Mar2018-v1_MINIAOD.root

#python ssbanalyzer_cfg.py print          \
cmsRun ssbanalyzer_cfg.py print          \
    label=DoubleMuon			 \
    outputFile=SSBTree.root    		 \
    globalTag=94X_dataRun2_v6		 \
    maxEvents=${EVENTS}        		 \
    inputFiles=${MYFILE}

