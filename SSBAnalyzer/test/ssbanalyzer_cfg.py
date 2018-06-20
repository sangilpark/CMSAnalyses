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
isSys=False


#data
print " >> label:: ", label
if label in [ 'SingleElectron', 'DoubleEG', 'SingleMuon', 'DoubleMuon', 'MuonEG', 'MET', 'SinglePhoton']:
    print " >> DATA:: ", label
    isMC = False


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
    print ("Running on Data ...globaltag : "+str(process.GlobalTag.globaltag))
else:
    #fname = 'file:/u/user/sangilpark/WorkDir/CP_Violation/TestSample/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_MINIAODSIM.root'
    #fname = 'file:/u/user/sangilpark/WorkDir/CP_Violation/TestSample/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1_MINIAODSIM.root'
    print ("Running on MC ...globaltag : "+str(process.GlobalTag.globaltag))


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

process.load("RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi")
process.load("RecoEgamma.PhotonIdentification.PhotonMVAValueMapProducer_cfi")

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce

my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
		 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
		 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',]


#add them to the VID producer
### Electron ID
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
                                                                
my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V1_TrueVtx_cff',
		    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff',]

#add them to the VID producer
### Photon ID
for idmod in my_phoid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

### Electron Cut for electron ID ###
process.selectedElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt>5 && abs(eta)")
)
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')
process.electronIDValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
process.electronRegressionValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')

### Photon Cut for photon ID ###
process.selectedPhotons = cms.EDFilter('PATPhotonSelector',
    src = cms.InputTag('slimmedPhotons'),
    cut = cms.string('pt>5 && abs(eta)')
)

process.egmPhotonIDs.physicsObjectSrc = cms.InputTag('selectedPhotons')
process.egmPhotonIsolation.srcToIsolate = cms.InputTag('selectedPhotons')
process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag('selectedPhotons')
process.photonRegressionValueMapProducer.srcMiniAOD = cms.InputTag('selectedPhotons')
process.photonMVAValueMapProducer.srcMiniAOD = cms.InputTag('selectedPhotons')


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



####################
### SSB Analyzer ###
####################

process.ssbanalyzer = cms.EDAnalyzer('SSBAnalyzer',
                                    bitsPat         = cms.InputTag("TriggerResults","","PAT"),
                                    isMCTag         = cms.bool(isMC),
                                    PDFInfoTag      = cms.InputTag("generator",""),
                                    FixPOWHEG	    = cms.untracked.string("NNPDF30_nlo_as_0118.LHgrid"),
                                    PDFSetNames     = cms.vstring('NNPDF30_nlo_as_0118.LHgrid','CT10nnlo.LHgrid'),#"cteq66.LHgrid"), #'MRST2006nnlo.LHgrid'),# 'NNPDF10_100.LHgrid'), 
                                    PDFCent         = cms.bool(True),
                                    PDFSys          = cms.bool(True),
                                    pvTag           = cms.InputTag("offlineSlimmedPrimaryVertices",""),
                                    genEvnTag	    = cms.InputTag("generator"),
                                    genLHETag	    = cms.InputTag("externalLHEProducer"),
                                    genParTag       = cms.InputTag("prunedGenParticles"),
                                    genJetTag       = cms.InputTag("slimmedGenJets","","PAT"),
                                    genJetReclusTag = cms.InputTag("slimmedGenJets","","PAT"), #FIXME
                                    #genBHadPlusMothersTag   = cms.InputTag("prunedGenParticles",), #FIXME
                                    #genBHadIndexTag   = cms.InputTag("slimmedGenJets","","PAT"), #FIXME
                                    #genBHadFlavourTag   = cms.InputTag("slimmedGenJets","","PAT"), #FIXME
                                    #genBHadFromTopWeakDecayTag   = cms.InputTag("slimmedGenJets","","PAT"), #FIXME
                                    #genBHadJetIndexTag   = cms.InputTag("slimmedGenJets","","PAT"), #FIXME

                                    genMETTag = cms.InputTag("slimmedMETs","","PAT"), #FIXME
                                    isSignal        = cms.bool(True),
                                    RhoTag          = cms.InputTag("fixedGridRhoFastjetAll"),
                                    puTag           = cms.InputTag("slimmedAddPileupInfo",""),
                                    trigList        = cms.vstring(
                                                        ### SingleMuon ###
                                                        'HLT_IsoMu24_v',
                                                        ### SingleElectron ###
                                                        'HLT_Ele32_WPTight_Gsf_v',
                                                        'HLT_Ele35_WPTight_Gsf_v',
                                                        'HLT_Ele38_WPTight_Gsf_v',
                                                        ### DoubleMuon ###
                                                        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_',
                                                        ### DoubleElectron ###
                                                        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                                        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                                        ### MuonElectron ###
                                                        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
                                                        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                                        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
                                                        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                                        ### AllHadronic ### 
                                                        'HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v',
                                                        'HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v',
                                                        ),
                                    bits             = cms.InputTag("TriggerResults","","HLT"),
                                    prescales        = cms.InputTag("patTrigger"),

                                    ismuSysTag       = cms.bool(isSys), # For Systemtic Study ...

                                    muTag            = cms.InputTag("slimmedMuons",""),
                                    muEnUpTag        = cms.InputTag("slimmedMuons",""),	#FIXME
                                    muEnDownTag      = cms.InputTag("slimmedMuons",""), #FIXME
                                    eleTag           = cms.InputTag("selectedElectrons",""), #FIXME
                                    electronPATInput = cms.InputTag("selectedElectrons",""), #FIXME
                                    eleEnUpTag       = cms.InputTag("selectedElectrons",""), #FIXME
                                    eleEnDownTag     = cms.InputTag("selectedElectrons",""), #FIXME

                                    effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt"),
                                    
				    eleVetoIdMap    = cms.InputTag( "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto"   ),
				    eleLooseIdMap   = cms.InputTag( "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose"  ),
				    eleMediumIdMap  = cms.InputTag( "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"   ),
				    eleTightIdMap   = cms.InputTag( "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight"    ),

                                    # ID decisions (mvaid)

				    mva_Iso_eleMediumMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-Iso-V1-wp90"),
				    mva_Iso_eleTightMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-Iso-V1-wp80"),
				    #mva_Iso_eleHZZIDMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-Iso-V1-wpLoose"),
				    mva_NoIso_eleMediumMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
				    mva_NoIso_eleTightMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
				    #mva_NoIso_eleHZZIDMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose"),
                                    
				    #
                                    # ValueMaps with MVA results
                                    #
				    mvaIsoValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
                                    mvaIsoCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
				    mvaNoIsoValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
                                    mvaNoIsoCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                                    #mvaValuesHZZMap  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
                                    #mvaCategoriesHZZMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Categories"),
                                    
                                    phoTag          = cms.InputTag("selectedPhotons",""),
                                    photonPATInput  = cms.InputTag("selectedPhotons",""), #FIXME
				    phoLooseIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V1-loose"  ),
				    phoMediumIdMap  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V1-medium" ),
				    phoTightIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Fall17-94X-V1-tight"  ), 

				    pho_mva_NontrigTightIdWP80Map = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v1-wp80"),
				    pho_mva_NontrigTightIdWP90Map = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v1-wp80"),

                                    pho_mvaWP80ValuesMap          = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Values"), #FIXME
                                    pho_mvaWP90ValuesMap          = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Values"), #FIXME
                                    full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                    phoChargedIsolation       = cms.InputTag('photonIDValueMapProducer:phoChargedIsolation'),
                                    phoNeutralHadronIsolation = cms.InputTag('photonIDValueMapProducer:phoNeutralHadronIsolation'),
                                    phoPhotonIsolation        = cms.InputTag('photonIDValueMapProducer:phoPhotonIsolation'),
                                    phoWorstChargedIsolation  = cms.InputTag('photonIDValueMapProducer:phoWorstChargedIsolation'),


                                    #effAreaChHadFile  = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased.txt"),
                                    #effAreaNeuHadFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt"),
                                    #effAreaPhoFile    = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased.txt"),
                                    effAreaChHadFile  = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_TrueVtx.txt"),
                                    effAreaNeuHadFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_TrueVtx.txt"),
                                    effAreaPhoFile    = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt"),

                                    pfCands         = cms.InputTag("packedPFCandidates"),
                                    bstag           = cms.InputTag("offlineBeamSpot"),
                                    convertag       = cms.InputTag("reducedEgamma","reducedConversions"),
                                    
                                    isjtcutTag = cms.bool(False), #For Jet
                                    jtTag = cms.InputTag("updatedPatJetsUpdatedJEC",""), #For Jet
                                    jtpuppiTag = cms.InputTag("slimmedJets",""), #For Jet
                                    PayLoadName = cms.string('AK4PFchs'),
                                    jer_useCondDB = cms.untracked.bool(True),  ### if jer_useCondDB == False, use external jer txts below
				    #phiResolMCFile = cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV10/Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt'),
                                    #phiResolDataFile = cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV6/Spring16_25nsV6_DATA_PhiResolution_AK4PFchs.txt'),
                                    #ptResolMCFile = cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV10/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),
                                    #ptResolDataFile = cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV6/Spring16_25nsV6_DATA_PtResolution_AK4PFchs.txt'),
                                    #ptResolSFFile = cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV10/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),
                                    
				    csvbjetTag = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                                    btagListTag        = cms.vstring(
                                                         'pfCombinedInclusiveSecondaryVertexV2BJetTags',
                                                         'softPFMuonBJetTags',
                                                         'softPFMuonByIP3dBJetTags',
                                                         'softPFElectronByPtBJetTags',
                                                         'softPFElectronBJetTags',
                                                         'softPFMuonByPtBJetTags',
                                                         'softPFElectronByIP3dBJetTags',
                                                         'softPFMuonByIP2dBJetTags',
                                                         'softPFElectronByIP2dBJetTags'
                                                      ),
                                    PFTightJetID = cms.PSet(
                                                       version = cms.string('RUNIISTARTUP'),
                                                       quality = cms.string('TIGHT')
                                                      ),
                                    PFLooseJetID = cms.PSet(
                                                       version = cms.string('RUNIISTARTUP'),
                                                       quality = cms.string('LOOSE')
                                                      ),
                                    metTag = cms.InputTag("slimmedMETs","","PAT"),
                                    metmucleancorTag = cms.InputTag("slimmedMETs","","PAT"),
)
process.p = cms.Path(
        process.selectedElectrons *
        process.selectedPhotons *
        process.egmGsfElectronIDSequence*
        process.egmPhotonIDSequence*
        process.photonIDValueMapProducer*
        process.electronIDValueMapProducer*
        process.photonMVAValueMapProducer *
        process.electronMVAValueMapProducer *
        process.JEC*
#        process.bfragStudy*
        process.ssbanalyzer
)
