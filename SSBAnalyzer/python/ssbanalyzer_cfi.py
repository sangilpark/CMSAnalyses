import FWCore.ParameterSet.Config as cms

####################
### SSB Analyzer ###
####################

ssbanalyzer = cms.EDAnalyzer('SSBAnalyzer',
	bitsPat         = cms.InputTag("TriggerResults","","PAT"),
	isMCTag         = cms.bool(True),
	PDFInfoTag      = cms.InputTag("generator",""),
	FixPOWHEG	= cms.untracked.string("NNPDF30_nlo_as_0118.LHgrid"),
	PDFSetNames     = cms.vstring('NNPDF30_nlo_as_0118.LHgrid','CT10nnlo.LHgrid'),#"cteq66.LHgrid"), #'MRST2006nnlo.LHgrid'),# 'NNPDF10_100.LHgrid'), 
	PDFCent         = cms.bool(True),
	PDFSys          = cms.bool(True),
	pvTag           = cms.InputTag("offlineSlimmedPrimaryVertices",""),
	genEvnTag	= cms.InputTag("generator"),
	genLHETag	= cms.InputTag("externalLHEProducer"),
	genParTag       = cms.InputTag("prunedGenParticles"),
	genJetTag       = cms.InputTag("slimmedGenJets","","PAT"),
	genMETTag	    = cms.InputTag("slimmedMETs","","PAT"), #FIXME

	isSignalTag     = cms.bool(False),
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

	ismuSysTag       = cms.bool(False), # For Systemtic Study ...

	muTag            = cms.InputTag("slimmedMuons",""),
	muEnUpTag        = cms.InputTag("slimmedMuons",""),	 #FIXME
	muEnDownTag      = cms.InputTag("slimmedMuons",""),	 #FIXME
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
	mva_NoIso_eleMediumMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
	mva_NoIso_eleTightMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
	
	# ValueMaps with MVA results
	mvaIsoValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
	mvaIsoCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
	mvaNoIsoValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
	mvaNoIsoCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
	
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
	jer_useCondDB = cms.untracked.bool(True),  ### Jet Energy Resolution, if jer_useCondDB == False, use external jer txts below
	#phiResolMCFile	= cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV10/Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt'),
	#phiResolDataFile	= cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV6/Spring16_25nsV6_DATA_PhiResolution_AK4PFchs.txt'),
	#ptResolMCFile	= cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV10/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),
	#ptResolDataFile	= cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV6/Spring16_25nsV6_DATA_PtResolution_AK4PFchs.txt'),
	#ptResolSFFile	= cms.FileInPath('CMSAnalyses/SSBAnalyzer/data/Spring16_25nsV10/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),
	
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
)
