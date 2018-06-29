import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import re
import sys

options = VarParsing.VarParsing('analysis')
options.register ('label',
                  'XXX',
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  'Label')

options.register ('globalTag',
                  '94X_dataRun2_v6',
                   VarParsing.VarParsing.multiplicity.singleton,
                   VarParsing.VarParsing.varType.string,
                  'GlobalTag')

#configurable options =======================================================================
#-------------------------------------------------------------------------------
# defaults
options.outputFile = 'SSBTree.root'
options.maxEvents = -1 # all events
#-------------------------------------------------------------------------------

options.parseArguments()


#===================================================================
# Define the CMSSW process
process = cms.Process("SSB")

label = options.label
globalTag = options.globalTag
isMC=True  #data/MC switch
isSignal=False 
isSys=False

#data, signal MC, background MC
print " >> label:: ", label
if label in [ 'SingleElectron', 'DoubleEG', 'SingleMuon', 'DoubleMuon', 'MuonEG', 'MET', 'SinglePhoton']:
    isMC = False
    print " >> DATA : ", label, " isMC : ", isMC
elif label in [ 'TTTo2L2Nu']:
    isSignal = True
    print " >> isSignal : " , isSignal , " isMC : ", isMC
else:
    print " >> isSignal : " , isSignal , " isMC : ", isMC


# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
   input = cms.untracked.int32(options.maxEvents),
)
########################
### Output filenames ###
########################
process.TFileService=cms.Service("TFileService",
        fileName=cms.string(options.outputFile),
        closeFileFast = cms.untracked.bool(True)
)

from Configuration.AlCa.GlobalTag import GlobalTag
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = globalTag

if not isMC: # Data
    #fname = 'file:/u/user/sangilpark/WorkDir/CP_Violation/TestSample/DoubleMuon_Run2017B-31Mar2018-v1_MINIAOD.root'
    print ("Running on Data ...globaltag : "+str(globalTag))
else:
    #fname = 'file:/u/user/sangilpark/WorkDir/CP_Violation/TestSample/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_MINIAODSIM.root'
    #fname = 'file:/u/user/sangilpark/WorkDir/CP_Violation/TestSample/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1_MINIAODSIM.root'
    print ("Running on MC ...globaltag : "+str(globalTag))


# Define the input source
process.source = cms.Source("PoolSource", 
    #fileNames = cms.untracked.vstring([ fname ])
    fileNames = cms.untracked.vstring(options.inputFiles)
)


#############################
### Random Number Service ###
#############################
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                  ssbanalyzer    = cms.PSet( initialSeed = cms.untracked.uint32(8675389),
                                             engineName = cms.untracked.string('TRandom3'),
                                                      ),
                                                   )

################################
### Electron & Photon -- ID ####
################################
### Applying the ID ######################################
process.load("RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi")
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

#process.load("RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi")
#process.load("RecoEgamma.PhotonIdentification.PhotonMVAValueMapProducer_cfi")

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
#switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
		 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
		 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',]

#my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V1_TrueVtx_cff',
#		    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff',]

#add them to the VID producer
### Electron ID
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

##add them to the VID producer
#### Photon ID
#for idmod in my_phoid_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

### Electron Cut for electron ID ###
process.selectedElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt>5 && abs(eta)")
)

#### Photon Cut for photon ID ###
#process.selectedPhotons = cms.EDFilter('PATPhotonSelector',
#    src = cms.InputTag('slimmedPhotons'),
#    cut = cms.string('pt>5 && abs(eta)')
#)

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')
process.electronIDValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
process.electronRegressionValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')

#process.egmPhotonIDs.physicsObjectSrc = cms.InputTag('selectedPhotons')
#process.egmPhotonIsolation.srcToIsolate = cms.InputTag('selectedPhotons')
#process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag('selectedPhotons')
#process.photonRegressionValueMapProducer.srcMiniAOD = cms.InputTag('selectedPhotons')
#process.photonMVAValueMapProducer.srcMiniAOD = cms.InputTag('selectedPhotons')

#############################
### Jet Energy Correction ###
#############################
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
)

process.JEC = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

"""
#############################################
### Jet Energy Resolution and SF using db ###
#############################################
process.load('Configuration.StandardSequences.Services_cff')
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *

process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        toGet = cms.VPSet(
            # Resolution
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Summer15_25nsV6_MC_PtResolution_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs_pt')
                ),

            # Scale factors
            cms.PSet(
                record = cms.string('JetResolutionScaleFactorRcd'),
                tag    = cms.string('JR_Summer15_25nsV6_MC_SF_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs')
                ),
            ),
        connect = cms.string('sqlite:Summer15_25nsV6.db')
        )

process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')
"""

## Following lines are for default MET for Type1 corrections.
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

process.load("CMSAnalyses.SSBAnalyzer.ssbanalyzer_cfi")
process.ssbanalyzer.isMCTag = isMC
process.ssbanalyzer.isSignalTag = isSignal

process.p = cms.Path(
        process.selectedElectrons *
#        process.selectedPhotons *
        process.egmGsfElectronIDSequence*
#        process.egmPhotonIDSequence*
#        process.photonIDValueMapProducer*
        process.electronIDValueMapProducer*
#        process.photonMVAValueMapProducer *
        process.electronMVAValueMapProducer *
        process.JEC*
#        process.bfragStudy*
        process.ssbanalyzer
)
