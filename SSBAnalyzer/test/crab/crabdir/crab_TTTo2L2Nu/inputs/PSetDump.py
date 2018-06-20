import FWCore.ParameterSet.Config as cms

process = cms.Process("SSB")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)
process.HFRecalParameterBlock = cms.PSet(
    HFdepthOneParameterA = cms.vdouble(0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
        0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
        0.058939, 0.125497),
    HFdepthOneParameterB = cms.vdouble(-4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
        2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
        0.000425, 0.000209),
    HFdepthTwoParameterA = cms.vdouble(0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
        0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
        0.051579, 0.086593),
    HFdepthTwoParameterB = cms.vdouble(-2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
        1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
        0.000157, -3e-06)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.mvaEleID_Fall17_iso_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
        -inf, -inf, -inf, -inf, -inf, 
        -1.0, -inf, -inf, -inf, -inf, 
        -inf, -inf, -0.06, -0.6, -0.2, 
        -inf, -inf, -inf, -inf, -inf),
    clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
        inf, inf, inf, 10.0, 200.0, 
        inf, inf, inf, inf, 20.0, 
        20.0, inf, 0.06, 0.6, 0.2, 
        inf, inf, inf, inf, inf),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    ebSplit = cms.double(0.8),
    ebeeSplit = cms.double(1.479),
    mvaName = cms.string('ElectronMVAEstimatorRun2Fall17Iso'),
    mvaTag = cms.string('V1'),
    ptSplit = cms.double(10.0),
    varNames = cms.vstring('ele_oldsigmaietaieta', 
        'ele_oldsigmaiphiiphi', 
        'ele_oldcircularity', 
        'ele_oldr9', 
        'ele_scletawidth', 
        'ele_sclphiwidth', 
        'ele_oldhe', 
        'ele_kfhits', 
        'ele_kfchi2', 
        'ele_gsfchi2', 
        'ele_fbrem', 
        'ele_gsfhits', 
        'ele_expected_inner_hits', 
        'ele_conversionVertexFitProbability', 
        'ele_ep', 
        'ele_eelepout', 
        'ele_IoEmIop', 
        'ele_deltaetain', 
        'ele_deltaphiin', 
        'ele_deltaetaseed', 
        'ele_pfPhotonIso', 
        'ele_pfChargedHadIso', 
        'ele_pfNeutralHadIso', 
        'rho', 
        'ele_psEoverEraw'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml')
)

process.mvaEleID_Fall17_noIso_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
        -inf, -inf, -inf, -inf, -inf, 
        -1.0, -inf, -inf, -inf, -inf, 
        -inf, -inf, -0.06, -0.6, -0.2, 
        -inf, -inf),
    clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
        inf, inf, inf, 10.0, 200.0, 
        inf, inf, inf, inf, 20.0, 
        20.0, inf, 0.06, 0.6, 0.2, 
        inf, inf),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    ebSplit = cms.double(0.8),
    ebeeSplit = cms.double(1.479),
    mvaName = cms.string('ElectronMVAEstimatorRun2Fall17NoIso'),
    mvaTag = cms.string('V1'),
    ptSplit = cms.double(10.0),
    varNames = cms.vstring('ele_oldsigmaietaieta', 
        'ele_oldsigmaiphiiphi', 
        'ele_oldcircularity', 
        'ele_oldr9', 
        'ele_scletawidth', 
        'ele_sclphiwidth', 
        'ele_oldhe', 
        'ele_kfhits', 
        'ele_kfchi2', 
        'ele_gsfchi2', 
        'ele_fbrem', 
        'ele_gsfhits', 
        'ele_expected_inner_hits', 
        'ele_conversionVertexFitProbability', 
        'ele_ep', 
        'ele_eelepout', 
        'ele_IoEmIop', 
        'ele_deltaetain', 
        'ele_deltaphiin', 
        'ele_deltaetaseed', 
        'rho', 
        'ele_psEoverEraw'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_BDT.weights.xml')
)

process.mvaEleID_PHYS14_PU20bx25_nonTrig_V1_producer_config = cms.PSet(
    mvaName = cms.string('ElectronMVAEstimatorRun2Phys14NonTrig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml')
)

process.mvaEleID_PHYS14_PU20bx25_nonTrig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.253, 0.081, -0.081, 0.965, 0.917, 
            0.683),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(False)
)

process.mvaEleID_PHYS14_PU20bx25_nonTrig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.483, -0.267, -0.323, 0.933, 0.825, 
            0.337),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(False)
)

process.mvaEleID_Spring15_25ns_Trig_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
)

process.mvaEleID_Spring15_25ns_Trig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
        mvaCuts = cms.vdouble(0.988153, 0.96791, 0.841729),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-Trig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_Trig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
        mvaCuts = cms.vdouble(0.972153, 0.922126, 0.610764),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-Trig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring15NonTrig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml')
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(0.287435, 0.221846, -0.303263, 0.967083, 0.929117, 
            0.726311),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.083313, -0.235222, -0.67099, 0.913286, 0.805013, 
            0.358969),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_wpLoose = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.265, -0.556, -0.551, -0.072, -0.286, 
            -0.267),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_50ns_Trig_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
    mvaTag = cms.string('50nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
)

process.mvaEleID_Spring15_50ns_Trig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Categories"),
        mvaCuts = cms.vdouble(0.981841, 0.946762, 0.79704),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-50ns-Trig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_50ns_Trig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Categories"),
        mvaCuts = cms.vdouble(0.953843, 0.849994, 0.514118),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-50ns-Trig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring16_GeneralPurpose_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring16GeneralPurpose'),
    mvaTag = cms.string('V1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml')
)

process.mvaEleID_Spring16_HZZ_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring16HZZ'),
    mvaTag = cms.string('V1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml')
)

process.mvaPhoID_RunIIFall17_v1_producer_config = cms.PSet(
    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
    mvaName = cms.string('PhotonMVAEstimatorRunIIFall17'),
    mvaTag = cms.string('v1'),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoIsoPtScalingCoeff = cms.vdouble(0.0035, 0.004),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationWithPVConstraint"),
    rho = cms.InputTag("fixedGridRhoAll"),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_barrel_BDT.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_endcap_BDT.weights.xml')
)

process.mvaPhoID_RunIIFall17_v1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1Categories"),
        mvaCuts = cms.vdouble(0.67, 0.54),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-RunIIFall17-v1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_RunIIFall17_v1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1Categories"),
        mvaCuts = cms.vdouble(0.27, 0.14),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-RunIIFall17-v1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_RunIIFall17_v1p1_producer_config = cms.PSet(
    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
    mvaName = cms.string('PhotonMVAEstimatorRunIIFall17'),
    mvaTag = cms.string('v1p1'),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoIsoPtScalingCoeff = cms.vdouble(0.0035, 0.004),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
    rho = cms.InputTag("fixedGridRhoAll"),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_barrel_BDT.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_endcap_BDT.weights.xml')
)

process.mvaPhoID_RunIIFall17_v1p1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1p1Categories"),
        mvaCuts = cms.vdouble(0.67, 0.54),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1p1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-RunIIFall17-v1p1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_RunIIFall17_v1p1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1p1Categories"),
        mvaCuts = cms.vdouble(0.27, 0.14),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1p1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-RunIIFall17-v1p1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_Spring15_25ns_nonTrig_V2p1_producer_config = cms.PSet(
    esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
    full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
    full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
    full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
    full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
    full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
    mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
    mvaTag = cms.string('25nsV2p1'),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    useValueMaps = cms.bool(False),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EB_V2.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EE_V2.weights.xml')
)

process.mvaPhoID_Spring15_25ns_nonTrig_V2p1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Categories"),
        mvaCuts = cms.vdouble(0.374, 0.336),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-Spring15-25ns-nonTrig-V2p1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_Spring15_50ns_nonTrig_V2p1_producer_config = cms.PSet(
    esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
    full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
    full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
    full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
    full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
    full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
    mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
    mvaTag = cms.string('50nsV2p1'),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    useValueMaps = cms.bool(False),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EB_V2.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EE_V2.weights.xml')
)

process.mvaPhoID_Spring15_50ns_nonTrig_V2p1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring15NonTrig50nsV2p1Categories"),
        mvaCuts = cms.vdouble(0.29538, 0.45837),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring15NonTrig50nsV2p1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-Spring15-50ns-nonTrig-V2p1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_Spring16_nonTrig_V1_producer_config = cms.PSet(
    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased_3bins.txt'),
    mvaName = cms.string('PhotonMVAEstimatorRun2Spring16NonTrig'),
    mvaTag = cms.string('V1'),
    phoChargedIsolation = cms.InputTag("egmPhotonIsolation","h+-DR030-"),
    phoIsoCutoff = cms.double(2.5),
    phoIsoPtScalingCoeff = cms.vdouble(0.0053, 0.0034),
    phoPhotonIsolation = cms.InputTag("egmPhotonIsolation","gamma-DR030-"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationWithConeVeto"),
    rho = cms.InputTag("fixedGridRhoAll"),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EB_V3.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EE_V3.weights.xml')
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.mvaConfigsForEleProducer = cms.VPSet(cms.PSet(
    mvaName = cms.string('ElectronMVAEstimatorRun2Phys14NonTrig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml')
), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring15NonTrig'),
        mvaTag = cms.string('25nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
        mvaTag = cms.string('50nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
        mvaTag = cms.string('25nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring16HZZ'),
        mvaTag = cms.string('V1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring16GeneralPurpose'),
        mvaTag = cms.string('V1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
            -inf, -inf, -inf, -inf, -inf, 
            -1.0, -inf, -inf, -inf, -inf, 
            -inf, -inf, -0.06, -0.6, -0.2, 
            -inf, -inf),
        clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
            inf, inf, inf, 10.0, 200.0, 
            inf, inf, inf, inf, 20.0, 
            20.0, inf, 0.06, 0.6, 0.2, 
            inf, inf),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        ebSplit = cms.double(0.8),
        ebeeSplit = cms.double(1.479),
        mvaName = cms.string('ElectronMVAEstimatorRun2Fall17NoIso'),
        mvaTag = cms.string('V1'),
        ptSplit = cms.double(10.0),
        varNames = cms.vstring('ele_oldsigmaietaieta', 
            'ele_oldsigmaiphiiphi', 
            'ele_oldcircularity', 
            'ele_oldr9', 
            'ele_scletawidth', 
            'ele_sclphiwidth', 
            'ele_oldhe', 
            'ele_kfhits', 
            'ele_kfchi2', 
            'ele_gsfchi2', 
            'ele_fbrem', 
            'ele_gsfhits', 
            'ele_expected_inner_hits', 
            'ele_conversionVertexFitProbability', 
            'ele_ep', 
            'ele_eelepout', 
            'ele_IoEmIop', 
            'ele_deltaetain', 
            'ele_deltaphiin', 
            'ele_deltaetaseed', 
            'rho', 
            'ele_psEoverEraw'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
            -inf, -inf, -inf, -inf, -inf, 
            -1.0, -inf, -inf, -inf, -inf, 
            -inf, -inf, -0.06, -0.6, -0.2, 
            -inf, -inf, -inf, -inf, -inf),
        clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
            inf, inf, inf, 10.0, 200.0, 
            inf, inf, inf, inf, 20.0, 
            20.0, inf, 0.06, 0.6, 0.2, 
            inf, inf, inf, inf, inf),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        ebSplit = cms.double(0.8),
        ebeeSplit = cms.double(1.479),
        mvaName = cms.string('ElectronMVAEstimatorRun2Fall17Iso'),
        mvaTag = cms.string('V1'),
        ptSplit = cms.double(10.0),
        varNames = cms.vstring('ele_oldsigmaietaieta', 
            'ele_oldsigmaiphiiphi', 
            'ele_oldcircularity', 
            'ele_oldr9', 
            'ele_scletawidth', 
            'ele_sclphiwidth', 
            'ele_oldhe', 
            'ele_kfhits', 
            'ele_kfchi2', 
            'ele_gsfchi2', 
            'ele_fbrem', 
            'ele_gsfhits', 
            'ele_expected_inner_hits', 
            'ele_conversionVertexFitProbability', 
            'ele_ep', 
            'ele_eelepout', 
            'ele_IoEmIop', 
            'ele_deltaetain', 
            'ele_deltaphiin', 
            'ele_deltaetaseed', 
            'ele_pfPhotonIso', 
            'ele_pfChargedHadIso', 
            'ele_pfNeutralHadIso', 
            'rho', 
            'ele_psEoverEraw'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml')
    ))

process.mvaConfigsForPhoProducer = cms.VPSet(cms.PSet(
    esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
    full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
    full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
    full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
    full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
    full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
    mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
    mvaTag = cms.string('50nsV2p1'),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    useValueMaps = cms.bool(False),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EB_V2.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EE_V2.weights.xml')
), 
    cms.PSet(
        esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
        full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
        full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
        full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
        full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
        full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
        full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
        mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
        mvaTag = cms.string('25nsV2p1'),
        phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
        phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
        phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        useValueMaps = cms.bool(False),
        weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EB_V2.weights.xml', 
            'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EE_V2.weights.xml')
    ), 
    cms.PSet(
        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased_3bins.txt'),
        mvaName = cms.string('PhotonMVAEstimatorRun2Spring16NonTrig'),
        mvaTag = cms.string('V1'),
        phoChargedIsolation = cms.InputTag("egmPhotonIsolation","h+-DR030-"),
        phoIsoCutoff = cms.double(2.5),
        phoIsoPtScalingCoeff = cms.vdouble(0.0053, 0.0034),
        phoPhotonIsolation = cms.InputTag("egmPhotonIsolation","gamma-DR030-"),
        phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationWithConeVeto"),
        rho = cms.InputTag("fixedGridRhoAll"),
        weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EB_V3.weights.xml', 
            'RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EE_V3.weights.xml')
    ), 
    cms.PSet(
        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
        mvaName = cms.string('PhotonMVAEstimatorRunIIFall17'),
        mvaTag = cms.string('v1'),
        phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
        phoIsoPtScalingCoeff = cms.vdouble(0.0035, 0.004),
        phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
        phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationWithPVConstraint"),
        rho = cms.InputTag("fixedGridRhoAll"),
        weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_barrel_BDT.weights.xml', 
            'RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_endcap_BDT.weights.xml')
    ), 
    cms.PSet(
        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
        mvaName = cms.string('PhotonMVAEstimatorRunIIFall17'),
        mvaTag = cms.string('v1p1'),
        phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
        phoIsoPtScalingCoeff = cms.vdouble(0.0035, 0.004),
        phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
        phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
        rho = cms.InputTag("fixedGridRhoAll"),
        weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_barrel_BDT.weights.xml', 
            'RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_endcap_BDT.weights.xml')
    ))

process.egmGsfElectronIDs = cms.EDProducer("VersionedGsfElectronIdProducer",
    physicsObjectIDs = cms.VPSet(cms.PSet(
        idDefinition = cms.PSet(
            cutFlow = cms.VPSet(cms.PSet(
                cutName = cms.string('MinPtCut'),
                isIgnored = cms.bool(False),
                minPt = cms.double(5.0),
                needsAdditionalProducts = cms.bool(False)
            ), 
                cms.PSet(
                    allowedEtaRanges = cms.VPSet(cms.PSet(
                        maxEta = cms.double(1.479),
                        minEta = cms.double(0.0)
                    ), 
                        cms.PSet(
                            maxEta = cms.double(2.5),
                            minEta = cms.double(1.479)
                        )),
                    cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False),
                    useAbsEta = cms.bool(True)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleDEtaInSeedCut'),
                    dEtaInSeedCutValueEB = cms.double(0.00387),
                    dEtaInSeedCutValueEE = cms.double(0.0072),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleDPhiInCut'),
                    dPhiInCutValueEB = cms.double(0.0716),
                    dPhiInCutValueEE = cms.double(0.147),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut'),
                    full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.0105),
                    full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0356),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelC0 = cms.double(0.05),
                    barrelCE = cms.double(1.12),
                    barrelCr = cms.double(0.0368),
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleHadronicOverEMEnergyScaledCut'),
                    endcapC0 = cms.double(0.0414),
                    endcapCE = cms.double(0.5),
                    endcapCr = cms.double(0.201),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    rho = cms.InputTag("fixedGridRhoFastjetAll")
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                    eInverseMinusPInverseCutValueEB = cms.double(0.129),
                    eInverseMinusPInverseCutValueEE = cms.double(0.0875),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleEffAreaPFIsoCut'),
                    effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt'),
                    isIgnored = cms.bool(False),
                    isRelativeIso = cms.bool(True),
                    isoCutEBHighPt = cms.double(0.133),
                    isoCutEBLowPt = cms.double(0.133),
                    isoCutEEHighPt = cms.double(0.146),
                    isoCutEELowPt = cms.double(0.146),
                    needsAdditionalProducts = cms.bool(True),
                    ptCutOff = cms.double(20.0),
                    rho = cms.InputTag("fixedGridRhoFastjetAll")
                ), 
                cms.PSet(
                    beamspotSrc = cms.InputTag("offlineBeamSpot"),
                    conversionSrc = cms.InputTag("allConversions"),
                    conversionSrcMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
                    cutName = cms.string('GsfEleConversionVetoCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleMissingHitsCut'),
                    isIgnored = cms.bool(False),
                    maxMissingHitsEB = cms.uint32(1),
                    maxMissingHitsEE = cms.uint32(1),
                    needsAdditionalProducts = cms.bool(False)
                )),
            idName = cms.string('cutBasedElectronID-Fall17-94X-V1-loose'),
            isPOGApproved = cms.untracked.bool(True)
        ),
        idMD5 = cms.string('0b8456d622494441fe713a6858e0f7c1'),
        isPOGApproved = cms.untracked.bool(True)
    ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDEtaInSeedCut'),
                        dEtaInSeedCutValueEB = cms.double(0.00365),
                        dEtaInSeedCutValueEE = cms.double(0.00625),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDPhiInCut'),
                        dPhiInCutValueEB = cms.double(0.0588),
                        dPhiInCutValueEE = cms.double(0.0355),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut'),
                        full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.0105),
                        full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0309),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelC0 = cms.double(0.026),
                        barrelCE = cms.double(1.12),
                        barrelCr = cms.double(0.0368),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleHadronicOverEMEnergyScaledCut'),
                        endcapC0 = cms.double(0.026),
                        endcapCE = cms.double(0.5),
                        endcapCr = cms.double(0.201),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll")
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                        eInverseMinusPInverseCutValueEB = cms.double(0.0327),
                        eInverseMinusPInverseCutValueEE = cms.double(0.0335),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEffAreaPFIsoCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt'),
                        isIgnored = cms.bool(False),
                        isRelativeIso = cms.bool(True),
                        isoCutEBHighPt = cms.double(0.0718),
                        isoCutEBLowPt = cms.double(0.0718),
                        isoCutEEHighPt = cms.double(0.143),
                        isoCutEELowPt = cms.double(0.143),
                        needsAdditionalProducts = cms.bool(True),
                        ptCutOff = cms.double(20.0),
                        rho = cms.InputTag("fixedGridRhoFastjetAll")
                    ), 
                    cms.PSet(
                        beamspotSrc = cms.InputTag("offlineBeamSpot"),
                        conversionSrc = cms.InputTag("allConversions"),
                        conversionSrcMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
                        cutName = cms.string('GsfEleConversionVetoCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleMissingHitsCut'),
                        isIgnored = cms.bool(False),
                        maxMissingHitsEB = cms.uint32(1),
                        maxMissingHitsEE = cms.uint32(1),
                        needsAdditionalProducts = cms.bool(False)
                    )),
                idName = cms.string('cutBasedElectronID-Fall17-94X-V1-medium'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('a238ee70910de53d36866e89768500e9'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDEtaInSeedCut'),
                        dEtaInSeedCutValueEB = cms.double(0.00353),
                        dEtaInSeedCutValueEE = cms.double(0.00567),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDPhiInCut'),
                        dPhiInCutValueEB = cms.double(0.0499),
                        dPhiInCutValueEE = cms.double(0.0165),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut'),
                        full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.0104),
                        full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0305),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelC0 = cms.double(0.026),
                        barrelCE = cms.double(1.12),
                        barrelCr = cms.double(0.0368),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleHadronicOverEMEnergyScaledCut'),
                        endcapC0 = cms.double(0.026),
                        endcapCE = cms.double(0.5),
                        endcapCr = cms.double(0.201),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll")
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                        eInverseMinusPInverseCutValueEB = cms.double(0.0278),
                        eInverseMinusPInverseCutValueEE = cms.double(0.0158),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEffAreaPFIsoCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt'),
                        isIgnored = cms.bool(False),
                        isRelativeIso = cms.bool(True),
                        isoCutEBHighPt = cms.double(0.0361),
                        isoCutEBLowPt = cms.double(0.0361),
                        isoCutEEHighPt = cms.double(0.094),
                        isoCutEELowPt = cms.double(0.094),
                        needsAdditionalProducts = cms.bool(True),
                        ptCutOff = cms.double(20.0),
                        rho = cms.InputTag("fixedGridRhoFastjetAll")
                    ), 
                    cms.PSet(
                        beamspotSrc = cms.InputTag("offlineBeamSpot"),
                        conversionSrc = cms.InputTag("allConversions"),
                        conversionSrcMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
                        cutName = cms.string('GsfEleConversionVetoCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleMissingHitsCut'),
                        isIgnored = cms.bool(False),
                        maxMissingHitsEB = cms.uint32(1),
                        maxMissingHitsEE = cms.uint32(1),
                        needsAdditionalProducts = cms.bool(False)
                    )),
                idName = cms.string('cutBasedElectronID-Fall17-94X-V1-tight'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('4acb2d2796efde7fba75380ce8823fc2'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDEtaInSeedCut'),
                        dEtaInSeedCutValueEB = cms.double(0.00523),
                        dEtaInSeedCutValueEE = cms.double(0.00984),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDPhiInCut'),
                        dPhiInCutValueEB = cms.double(0.159),
                        dPhiInCutValueEE = cms.double(0.157),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut'),
                        full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.0128),
                        full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0445),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelC0 = cms.double(0.05),
                        barrelCE = cms.double(1.12),
                        barrelCr = cms.double(0.0368),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleHadronicOverEMEnergyScaledCut'),
                        endcapC0 = cms.double(0.05),
                        endcapCE = cms.double(0.5),
                        endcapCr = cms.double(0.201),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll")
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                        eInverseMinusPInverseCutValueEB = cms.double(0.193),
                        eInverseMinusPInverseCutValueEE = cms.double(0.0962),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEffAreaPFIsoCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt'),
                        isIgnored = cms.bool(False),
                        isRelativeIso = cms.bool(True),
                        isoCutEBHighPt = cms.double(0.168),
                        isoCutEBLowPt = cms.double(0.168),
                        isoCutEEHighPt = cms.double(0.185),
                        isoCutEELowPt = cms.double(0.185),
                        needsAdditionalProducts = cms.bool(True),
                        ptCutOff = cms.double(20.0),
                        rho = cms.InputTag("fixedGridRhoFastjetAll")
                    ), 
                    cms.PSet(
                        beamspotSrc = cms.InputTag("offlineBeamSpot"),
                        conversionSrc = cms.InputTag("allConversions"),
                        conversionSrcMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
                        cutName = cms.string('GsfEleConversionVetoCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleMissingHitsCut'),
                        isIgnored = cms.bool(False),
                        maxMissingHitsEB = cms.uint32(2),
                        maxMissingHitsEE = cms.uint32(3),
                        needsAdditionalProducts = cms.bool(False)
                    )),
                idName = cms.string('cutBasedElectronID-Fall17-94X-V1-veto'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('43be9b381a8d9b0910b7f81a5ad8ff3a'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVAExpoScalingCut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                    mvaCuts = cms.vdouble(0.953024095656, 2.7591425841, 0.466964471855, 0.933656476396, 2.70927628427, 
                        0.335122865992, 0.931313368837, 1.58219348007, 3.88894626197, 0.982526856494, 
                        8.70260145586, 1.19748615966, 0.972750945793, 8.17952563102, 1.71117550947, 
                        0.956261953954, 8.10984536628, 3.01392769913),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-noIso-V1-wp80'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVAExpoScalingCut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                    mvaCuts = cms.vdouble(0.916511282697, 2.73817035551, 1.03549199648, 0.865573832222, 2.40279446526, 
                        0.797561561328, -3016.03505523, -52140.6185633, -3016.30293872, 0.961654281613, 
                        8.75794383789, 3.13902003216, 0.931925801143, 8.84605743257, 3.59850637933, 
                        0.8899260781, 10.1242341159, 4.35279125072),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-noIso-V1-wp90'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                    mvaCuts = cms.vdouble(-0.132858672938, -0.317653009588, -0.0799205914719, -0.856871961305, -0.810764214158, 
                        -0.717926593302),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-noIso-V1-wpLoose'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVAExpoScalingCut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
                    mvaCuts = cms.vdouble(0.972550955975, 2.97659326151, 0.26538587364, 0.95080381416, 2.66335005587, 
                        0.235582049926, 0.93650371676, 1.57654423239, 3.06701528922, 0.989656208772, 
                        10.342490512, 0.402041564174, 0.981923265653, 9.05548836482, 0.772674931169, 
                        0.962509820174, 8.42589315557, 2.29161526151),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-iso-V1-wp80'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVAExpoScalingCut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
                    mvaCuts = cms.vdouble(0.93870703961, 2.65255852282, 0.822264716415, 0.894880292568, 2.76456703588, 
                        0.41233812187, -1830.85836611, -36578.1105538, -1831.20835781, 0.971767483761, 
                        8.9128509851, 1.97124149404, 0.945874502327, 8.83104420393, 2.40849932041, 
                        0.897911201209, 9.81408214417, 4.17158169489),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-iso-V1-wp90'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
                    mvaCuts = cms.vdouble(-0.0956408614642, -0.282299169819, -0.0546668229696, -0.833466688584, -0.767700024757, 
                        -0.691730599565),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-iso-V1-wpLoose'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        )),
    physicsObjectSrc = cms.InputTag("selectedElectrons")
)


process.egmPhotonIDs = cms.EDProducer("VersionedPhotonIdProducer",
    physicsObjectIDs = cms.VPSet(cms.PSet(
        idDefinition = cms.PSet(
            cutFlow = cms.VPSet(cms.PSet(
                cutName = cms.string('MinPtCut'),
                isIgnored = cms.bool(False),
                minPt = cms.double(5.0),
                needsAdditionalProducts = cms.bool(False)
            ), 
                cms.PSet(
                    allowedEtaRanges = cms.VPSet(cms.PSet(
                        maxEta = cms.double(1.479),
                        minEta = cms.double(0.0)
                    ), 
                        cms.PSet(
                            maxEta = cms.double(2.5),
                            minEta = cms.double(1.479)
                        )),
                    cutName = cms.string('PhoSCEtaMultiRangeCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False),
                    useAbsEta = cms.bool(True)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoSingleTowerHadOverEmCut'),
                    hadronicOverEMCutValueEB = cms.double(0.105),
                    hadronicOverEMCutValueEE = cms.double(0.029),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoFull5x5SigmaIEtaIEtaCut'),
                    cutValueEB = cms.double(0.0103),
                    cutValueEE = cms.double(0.0276),
                    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    C1_EB = cms.double(2.839),
                    C1_EE = cms.double(2.15),
                    C2_EB = cms.double(0),
                    C2_EE = cms.double(0.0),
                    anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoAnyPFIsoWithEACut'),
                    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_TrueVtx.txt'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                    useRelativeIso = cms.bool(False)
                ), 
                cms.PSet(
                    C1_EB = cms.double(9.188),
                    C1_EE = cms.double(10.471),
                    C2_EB = cms.double(0.0126),
                    C2_EE = cms.double(0.0119),
                    C3_EB = cms.double(2.6e-05),
                    C3_EE = cms.double(2.5e-05),
                    anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoAnyPFIsoWithEAAndQuadScalingCut'),
                    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_TrueVtx.txt'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                    useRelativeIso = cms.bool(False)
                ), 
                cms.PSet(
                    C1_EB = cms.double(2.956),
                    C1_EE = cms.double(4.895),
                    C2_EB = cms.double(0.0035),
                    C2_EE = cms.double(0.004),
                    anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoAnyPFIsoWithEACut'),
                    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                    useRelativeIso = cms.bool(False)
                )),
            idName = cms.string('cutBasedPhotonID-Fall17-94X-V1-loose'),
            isPOGApproved = cms.untracked.bool(True)
        ),
        idMD5 = cms.string('45515ee95e01fa36972ff7ba69186c97'),
        isPOGApproved = cms.untracked.bool(True)
    ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('PhoSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoSingleTowerHadOverEmCut'),
                        hadronicOverEMCutValueEB = cms.double(0.035),
                        hadronicOverEMCutValueEE = cms.double(0.027),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoFull5x5SigmaIEtaIEtaCut'),
                        cutValueEB = cms.double(0.0103),
                        cutValueEE = cms.double(0.0271),
                        full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(1.416),
                        C1_EE = cms.double(1.012),
                        C2_EB = cms.double(0.0),
                        C2_EE = cms.double(0.0),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEACut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_TrueVtx.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(2.491),
                        C1_EE = cms.double(9.131),
                        C2_EB = cms.double(0.0126),
                        C2_EE = cms.double(0.0119),
                        C3_EB = cms.double(2.6e-05),
                        C3_EE = cms.double(2.5e-05),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEAAndQuadScalingCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_TrueVtx.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(2.952),
                        C1_EE = cms.double(4.095),
                        C2_EB = cms.double(0.004),
                        C2_EE = cms.double(0.004),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEACut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    )),
                idName = cms.string('cutBasedPhotonID-Fall17-94X-V1-medium'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('772f7921fa146b630e4dbe79e475a421'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('PhoSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoSingleTowerHadOverEmCut'),
                        hadronicOverEMCutValueEB = cms.double(0.02),
                        hadronicOverEMCutValueEE = cms.double(0.025),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoFull5x5SigmaIEtaIEtaCut'),
                        cutValueEB = cms.double(0.0103),
                        cutValueEE = cms.double(0.0271),
                        full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(1.158),
                        C1_EE = cms.double(0.575),
                        C2_EB = cms.double(0.0),
                        C2_EE = cms.double(0.0),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEACut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_TrueVtx.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(1.267),
                        C1_EE = cms.double(8.916),
                        C2_EB = cms.double(0.0126),
                        C2_EE = cms.double(0.0119),
                        C3_EB = cms.double(2.6e-05),
                        C3_EE = cms.double(2.5e-05),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEAAndQuadScalingCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_TrueVtx.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(2.065),
                        C1_EE = cms.double(3.272),
                        C2_EB = cms.double(0.0035),
                        C2_EE = cms.double(0.004),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEACut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    )),
                idName = cms.string('cutBasedPhotonID-Fall17-94X-V1-tight'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('e260fee6f9011fb13ff56d45cccd21c5'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('PhoMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1Categories"),
                    mvaCuts = cms.vdouble(0.67, 0.54),
                    mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaPhoID-RunIIFall17-v1-wp80'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('7e48c47329d7d1eb889100ed03a02ba9'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('PhoMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1Categories"),
                    mvaCuts = cms.vdouble(0.27, 0.14),
                    mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaPhoID-RunIIFall17-v1-wp90'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('834dd792692b6a62786bd1caa6b53a68'),
            isPOGApproved = cms.untracked.bool(True)
        )),
    physicsObjectSrc = cms.InputTag("selectedPhotons")
)


process.egmPhotonIsolation = cms.EDProducer("CITKPFIsolationSumProducer",
    isolationConeDefinitions = cms.VPSet(cms.PSet(
        coneSize = cms.double(0.3),
        isolateAgainst = cms.string('h+'),
        isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
        miniAODVertexCodes = cms.vuint32(2, 3),
        particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
        vertexIndex = cms.int32(0)
    ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('h0'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
            vertexIndex = cms.int32(0)
        ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('gamma'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
            vertexIndex = cms.int32(0)
        )),
    srcForIsolationCone = cms.InputTag("packedPFCandidates"),
    srcToIsolate = cms.InputTag("selectedPhotons")
)


process.electronIDValueMapProducer = cms.EDProducer("ElectronIDValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("selectedElectrons")
)


process.electronMVAValueMapProducer = cms.EDProducer("ElectronMVAValueMapProducer",
    mvaConfigurations = cms.VPSet(cms.PSet(
        mvaName = cms.string('ElectronMVAEstimatorRun2Phys14NonTrig'),
        mvaTag = cms.string('25nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml')
    ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring15NonTrig'),
            mvaTag = cms.string('25nsV1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
            mvaTag = cms.string('50nsV1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
            mvaTag = cms.string('25nsV1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring16HZZ'),
            mvaTag = cms.string('V1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring16GeneralPurpose'),
            mvaTag = cms.string('V1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
                -inf, -inf, -inf, -inf, -inf, 
                -1.0, -inf, -inf, -inf, -inf, 
                -inf, -inf, -0.06, -0.6, -0.2, 
                -inf, -inf),
            clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
                inf, inf, inf, 10.0, 200.0, 
                inf, inf, inf, inf, 20.0, 
                20.0, inf, 0.06, 0.6, 0.2, 
                inf, inf),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            ebSplit = cms.double(0.8),
            ebeeSplit = cms.double(1.479),
            mvaName = cms.string('ElectronMVAEstimatorRun2Fall17NoIso'),
            mvaTag = cms.string('V1'),
            ptSplit = cms.double(10.0),
            varNames = cms.vstring('ele_oldsigmaietaieta', 
                'ele_oldsigmaiphiiphi', 
                'ele_oldcircularity', 
                'ele_oldr9', 
                'ele_scletawidth', 
                'ele_sclphiwidth', 
                'ele_oldhe', 
                'ele_kfhits', 
                'ele_kfchi2', 
                'ele_gsfchi2', 
                'ele_fbrem', 
                'ele_gsfhits', 
                'ele_expected_inner_hits', 
                'ele_conversionVertexFitProbability', 
                'ele_ep', 
                'ele_eelepout', 
                'ele_IoEmIop', 
                'ele_deltaetain', 
                'ele_deltaphiin', 
                'ele_deltaetaseed', 
                'rho', 
                'ele_psEoverEraw'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
                -inf, -inf, -inf, -inf, -inf, 
                -1.0, -inf, -inf, -inf, -inf, 
                -inf, -inf, -0.06, -0.6, -0.2, 
                -inf, -inf, -inf, -inf, -inf),
            clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
                inf, inf, inf, 10.0, 200.0, 
                inf, inf, inf, inf, 20.0, 
                20.0, inf, 0.06, 0.6, 0.2, 
                inf, inf, inf, inf, inf),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            ebSplit = cms.double(0.8),
            ebeeSplit = cms.double(1.479),
            mvaName = cms.string('ElectronMVAEstimatorRun2Fall17Iso'),
            mvaTag = cms.string('V1'),
            ptSplit = cms.double(10.0),
            varNames = cms.vstring('ele_oldsigmaietaieta', 
                'ele_oldsigmaiphiiphi', 
                'ele_oldcircularity', 
                'ele_oldr9', 
                'ele_scletawidth', 
                'ele_sclphiwidth', 
                'ele_oldhe', 
                'ele_kfhits', 
                'ele_kfchi2', 
                'ele_gsfchi2', 
                'ele_fbrem', 
                'ele_gsfhits', 
                'ele_expected_inner_hits', 
                'ele_conversionVertexFitProbability', 
                'ele_ep', 
                'ele_eelepout', 
                'ele_IoEmIop', 
                'ele_deltaetain', 
                'ele_deltaphiin', 
                'ele_deltaetaseed', 
                'ele_pfPhotonIso', 
                'ele_pfChargedHadIso', 
                'ele_pfNeutralHadIso', 
                'rho', 
                'ele_psEoverEraw'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml')
        )),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("selectedElectrons")
)


process.electronRegressionValueMapProducer = cms.EDProducer("ElectronRegressionValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("selectedElectrons"),
    useFull5x5 = cms.bool(False)
)


process.patJetCorrFactorsUpdatedJEC = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring('L1FastJet', 
        'L2Relative', 
        'L3Absolute', 
        'L2L3Residual'),
    payload = cms.string('AK4PFchs'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedJets"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.photonIDValueMapProducer = cms.EDProducer("PhotonIDValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
    pfCandidates = cms.InputTag("particleFlow"),
    pfCandidatesMiniAOD = cms.InputTag("packedPFCandidates"),
    src = cms.InputTag("gedPhotons"),
    srcMiniAOD = cms.InputTag("selectedPhotons"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    verticesMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.photonMVAValueMapProducer = cms.EDProducer("PhotonMVAValueMapProducer",
    mvaConfigurations = cms.VPSet(cms.PSet(
        esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
        full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
        full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
        full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
        full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
        full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
        full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
        mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
        mvaTag = cms.string('50nsV2p1'),
        phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
        phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
        phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        useValueMaps = cms.bool(False),
        weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EB_V2.weights.xml', 
            'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EE_V2.weights.xml')
    ), 
        cms.PSet(
            esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
            full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
            full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
            full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
            full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
            full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
            full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
            mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
            mvaTag = cms.string('25nsV2p1'),
            phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
            phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
            phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
            rho = cms.InputTag("fixedGridRhoFastjetAll"),
            useValueMaps = cms.bool(False),
            weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EB_V2.weights.xml', 
                'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EE_V2.weights.xml')
        ), 
        cms.PSet(
            effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased_3bins.txt'),
            mvaName = cms.string('PhotonMVAEstimatorRun2Spring16NonTrig'),
            mvaTag = cms.string('V1'),
            phoChargedIsolation = cms.InputTag("egmPhotonIsolation","h+-DR030-"),
            phoIsoCutoff = cms.double(2.5),
            phoIsoPtScalingCoeff = cms.vdouble(0.0053, 0.0034),
            phoPhotonIsolation = cms.InputTag("egmPhotonIsolation","gamma-DR030-"),
            phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationWithConeVeto"),
            rho = cms.InputTag("fixedGridRhoAll"),
            weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EB_V3.weights.xml', 
                'RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EE_V3.weights.xml')
        ), 
        cms.PSet(
            effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
            mvaName = cms.string('PhotonMVAEstimatorRunIIFall17'),
            mvaTag = cms.string('v1'),
            phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
            phoIsoPtScalingCoeff = cms.vdouble(0.0035, 0.004),
            phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
            phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationWithPVConstraint"),
            rho = cms.InputTag("fixedGridRhoAll"),
            weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_barrel_BDT.weights.xml', 
                'RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_endcap_BDT.weights.xml')
        ), 
        cms.PSet(
            effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
            mvaName = cms.string('PhotonMVAEstimatorRunIIFall17'),
            mvaTag = cms.string('v1p1'),
            phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
            phoIsoPtScalingCoeff = cms.vdouble(0.0035, 0.004),
            phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
            phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
            rho = cms.InputTag("fixedGridRhoAll"),
            weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_barrel_BDT.weights.xml', 
                'RecoEgamma/PhotonIdentification/data/Fall17/HggPhoId_92X_endcap_BDT.weights.xml')
        )),
    src = cms.InputTag("gedPhotons"),
    srcMiniAOD = cms.InputTag("selectedPhotons")
)


process.photonRegressionValueMapProducer = cms.EDProducer("PhotonRegressionValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    src = cms.InputTag("gedPhotons"),
    srcMiniAOD = cms.InputTag("selectedPhotons"),
    useFull5x5 = cms.bool(False)
)


process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")


process.updatedPatJetsUpdatedJEC = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(False),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsUpdatedJEC")),
    jetSource = cms.InputTag("slimmedJets"),
    printWarning = cms.bool(True),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.selectedElectrons = cms.EDFilter("PATElectronSelector",
    cut = cms.string('pt>5 && abs(eta)'),
    src = cms.InputTag("slimmedElectrons")
)


process.selectedPhotons = cms.EDFilter("PATPhotonSelector",
    cut = cms.string('pt>5 && abs(eta)'),
    src = cms.InputTag("slimmedPhotons")
)


process.selectedUpdatedPatJetsUpdatedJEC = cms.EDFilter("PATJetSelector",
    cut = cms.string(''),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("updatedPatJetsUpdatedJEC")
)


process.ssbanalyzer = cms.EDAnalyzer("SSBAnalyzer",
    FixPOWHEG = cms.untracked.string('NNPDF30_nlo_as_0118.LHgrid'),
    PDFCent = cms.bool(True),
    PDFInfoTag = cms.InputTag("generator"),
    PDFSetNames = cms.vstring('NNPDF30_nlo_as_0118.LHgrid', 
        'CT10nnlo.LHgrid'),
    PDFSys = cms.bool(True),
    PFLooseJetID = cms.PSet(
        quality = cms.string('LOOSE'),
        version = cms.string('RUNIISTARTUP')
    ),
    PFTightJetID = cms.PSet(
        quality = cms.string('TIGHT'),
        version = cms.string('RUNIISTARTUP')
    ),
    PayLoadName = cms.string('AK4PFchs'),
    RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    bitsPat = cms.InputTag("TriggerResults","","PAT"),
    bstag = cms.InputTag("offlineBeamSpot"),
    btagListTag = cms.vstring('pfCombinedInclusiveSecondaryVertexV2BJetTags', 
        'softPFMuonBJetTags', 
        'softPFMuonByIP3dBJetTags', 
        'softPFElectronByPtBJetTags', 
        'softPFElectronBJetTags', 
        'softPFMuonByPtBJetTags', 
        'softPFElectronByIP3dBJetTags', 
        'softPFMuonByIP2dBJetTags', 
        'softPFElectronByIP2dBJetTags'),
    convertag = cms.InputTag("reducedEgamma","reducedConversions"),
    csvbjetTag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    effAreaChHadFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_TrueVtx.txt'),
    effAreaNeuHadFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_TrueVtx.txt'),
    effAreaPhoFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt'),
    effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt'),
    eleEnDownTag = cms.InputTag("selectedElectrons"),
    eleEnUpTag = cms.InputTag("selectedElectrons"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Fall17-94X-V1-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Fall17-94X-V1-medium"),
    eleTag = cms.InputTag("selectedElectrons"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Fall17-94X-V1-tight"),
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Fall17-94X-V1-veto"),
    electronPATInput = cms.InputTag("selectedElectrons"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
    genEvnTag = cms.InputTag("generator"),
    genJetReclusTag = cms.InputTag("slimmedGenJets","","PAT"),
    genJetTag = cms.InputTag("slimmedGenJets","","PAT"),
    genLHETag = cms.InputTag("externalLHEProducer"),
    genMETTag = cms.InputTag("slimmedMETs","","PAT"),
    genParTag = cms.InputTag("prunedGenParticles"),
    isMCTag = cms.bool(True),
    isSignal = cms.bool(True),
    isjtcutTag = cms.bool(False),
    ismuSysTag = cms.bool(False),
    jer_useCondDB = cms.untracked.bool(True),
    jtTag = cms.InputTag("updatedPatJetsUpdatedJEC"),
    jtpuppiTag = cms.InputTag("slimmedJets"),
    metTag = cms.InputTag("slimmedMETs","","PAT"),
    metmucleancorTag = cms.InputTag("slimmedMETs","","PAT"),
    muEnDownTag = cms.InputTag("slimmedMuons"),
    muEnUpTag = cms.InputTag("slimmedMuons"),
    muTag = cms.InputTag("slimmedMuons"),
    mvaIsoCategoriesMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
    mvaIsoValuesMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Values"),
    mvaNoIsoCategoriesMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
    mvaNoIsoValuesMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
    mva_Iso_eleMediumMap = cms.InputTag("egmGsfElectronIDs","mvaEleID-Fall17-Iso-V1-wp90"),
    mva_Iso_eleTightMap = cms.InputTag("egmGsfElectronIDs","mvaEleID-Fall17-Iso-V1-wp80"),
    mva_NoIso_eleMediumMap = cms.InputTag("egmGsfElectronIDs","mvaEleID-Fall17-noIso-V1-wp90"),
    mva_NoIso_eleTightMap = cms.InputTag("egmGsfElectronIDs","mvaEleID-Fall17-noIso-V1-wp80"),
    pfCands = cms.InputTag("packedPFCandidates"),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoLooseIdMap = cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Fall17-94X-V1-loose"),
    phoMediumIdMap = cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Fall17-94X-V1-medium"),
    phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoTag = cms.InputTag("selectedPhotons"),
    phoTightIdMap = cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Fall17-94X-V1-tight"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
    pho_mvaWP80ValuesMap = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
    pho_mvaWP90ValuesMap = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
    pho_mva_NontrigTightIdWP80Map = cms.InputTag("egmPhotonIDs","mvaPhoID-RunIIFall17-v1-wp80"),
    pho_mva_NontrigTightIdWP90Map = cms.InputTag("egmPhotonIDs","mvaPhoID-RunIIFall17-v1-wp80"),
    photonPATInput = cms.InputTag("selectedPhotons"),
    prescales = cms.InputTag("patTrigger"),
    puTag = cms.InputTag("slimmedAddPileupInfo"),
    pvTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    trigList = cms.vstring('HLT_IsoMu24_v', 
        'HLT_Ele32_WPTight_Gsf_v', 
        'HLT_Ele35_WPTight_Gsf_v', 
        'HLT_Ele38_WPTight_Gsf_v', 
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_', 
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v', 
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v', 
        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v', 
        'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v', 
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v', 
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v', 
        'HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v', 
        'HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v')
)


process.DQMStore = cms.Service("DQMStore",
    LSbasedMode = cms.untracked.bool(False),
    collateHistograms = cms.untracked.bool(False),
    enableMultiThread = cms.untracked.bool(False),
    forceResetOnBeginLumi = cms.untracked.bool(False),
    referenceFileName = cms.untracked.string(''),
    verbose = cms.untracked.int32(0),
    verboseQT = cms.untracked.int32(0)
)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(100)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring('cout', 
        'cerr'),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    ssbanalyzer = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(8675389)
    )
)


process.TFileService = cms.Service("TFileService",
    closeFileFast = cms.untracked.bool(True),
    fileName = cms.string('SSBTree.root')
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring('HCAL', 
        'ZDC', 
        'CASTOR', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER')
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.MuonNumberingInitialization = cms.ESProducer("MuonNumberingInitialization")


process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(18268),
    version = cms.string('Parabolic')
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True),
    useDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(18268)
)


process.XMLFromDBSource = cms.ESProducer("XMLIdealGeometryESProducer",
    label = cms.string('Extended'),
    rootDDName = cms.string('cms:OCMS')
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(True),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiPixelQualityFromDbRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        ))
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(cms.PSet(
        Label = cms.untracked.string(''),
        NormalizationFactor = cms.untracked.double(1.0),
        Record = cms.string('SiStripApvGainRcd')
    ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiStripDetVOffRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('94X_mc2017_realistic_v14'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HepPDTESSource = cms.ESSource("HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/pythiaparticle.tbl')
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HBRecalibration = cms.bool(False),
    HBmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHB.txt'),
    HBreCalibCutoff = cms.double(20.0),
    HERecalibration = cms.bool(False),
    HEmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHE.txt'),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalParameterBlock = cms.PSet(
        HFdepthOneParameterA = cms.vdouble(0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
            0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
            0.058939, 0.125497),
        HFdepthOneParameterB = cms.vdouble(-4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
            2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
            0.000425, 0.000209),
        HFdepthTwoParameterA = cms.vdouble(0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
            0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
            0.051579, 0.086593),
        HFdepthTwoParameterB = cms.vdouble(-2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
            1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
            0.000157, -3e-06)
    ),
    HFRecalibration = cms.bool(False),
    SiPMCharacteristics = cms.VPSet(cms.PSet(
        crosstalk = cms.double(0.0),
        nonlin1 = cms.double(1.0),
        nonlin2 = cms.double(0.0),
        nonlin3 = cms.double(0.0),
        pixels = cms.int32(36000)
    ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(2500)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(0)
        )),
    hb = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.19),
        gainWidth = cms.vdouble(0.0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.285),
        pedestalWidth = cms.double(0.809),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.49, 1.8, 7.2, 37.9),
        qieSlope = cms.vdouble(0.912, 0.917, 0.922, 0.923),
        qieType = cms.int32(0),
        recoShape = cms.int32(105)
    ),
    hbUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.000439367311072),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(203),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(44.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.69e-11, 7.9e-11),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(150),
            intlumiToNeutrons = cms.double(367000000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(-5)
        ),
        recoShape = cms.int32(203)
    ),
    he = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.23),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.163),
        pedestalWidth = cms.double(0.9698),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.38, 2.0, 7.6, 39.6),
        qieSlope = cms.vdouble(0.912, 0.916, 0.92, 0.922),
        qieType = cms.int32(0),
        recoShape = cms.int32(105)
    ),
    heUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.000439367311072),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(203),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(44.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.69e-11, 7.9e-11),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(75),
            intlumiToNeutrons = cms.double(29200000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(5)
        ),
        recoShape = cms.int32(203)
    ),
    hf = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(9.354),
        pedestalWidth = cms.double(2.516),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(-0.87, 1.4, 7.8, -29.6),
        qieSlope = cms.vdouble(0.359, 0.358, 0.36, 0.367),
        qieType = cms.int32(0),
        recoShape = cms.int32(301)
    ),
    hfUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(13.33),
        pedestalWidth = cms.double(3.33),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(0.0697, -0.7405, 12.38, -671.9),
        qieSlope = cms.vdouble(0.297, 0.298, 0.298, 0.313),
        qieType = cms.int32(1),
        recoShape = cms.int32(301)
    ),
    ho = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.006, 0.0087),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(201),
        pedestal = cms.double(12.06),
        pedestalWidth = cms.double(0.6285),
        photoelectronsToAnalog = cms.double(4.0),
        qieOffset = cms.vdouble(-0.44, 1.4, 7.1, 38.5),
        qieSlope = cms.vdouble(0.907, 0.915, 0.92, 0.921),
        qieType = cms.int32(0),
        recoShape = cms.int32(201)
    ),
    iLumi = cms.double(-1.0),
    killHE = cms.bool(False),
    testHEPlan1 = cms.bool(False),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths'),
    useHBUpgrade = cms.bool(False),
    useHEUpgrade = cms.bool(False),
    useHFUpgrade = cms.bool(False),
    useHOUpgrade = cms.bool(True),
    useLayer0Weight = cms.bool(False)
)


process.prefer("es_hardcode")

process.patAlgosToolsTask = cms.Task(process.patJetCorrFactorsUpdatedJEC, process.selectedUpdatedPatJetsUpdatedJEC, process.updatedPatJetsUpdatedJEC)


process.egmPhotonIsolationMiniAODTask = cms.Task(process.egmPhotonIsolation)


process.egmGsfElectronIDTask = cms.Task(process.egmGsfElectronIDs, process.electronMVAValueMapProducer, process.electronRegressionValueMapProducer)


process.egmPhotonIDTask = cms.Task(process.egmPhotonIDs, process.egmPhotonIsolationMiniAODTask, process.photonIDValueMapProducer, process.photonMVAValueMapProducer, process.photonRegressionValueMapProducer)


process.JEC = cms.Sequence(process.patJetCorrFactorsUpdatedJEC+process.updatedPatJetsUpdatedJEC)


process.egmPhotonIDSequence = cms.Sequence(process.egmPhotonIDTask)


process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDTask)


process.egmPhotonIsolationMiniAODSequence = cms.Sequence(process.egmPhotonIsolationMiniAODTask)


process.p = cms.Path(process.selectedElectrons+process.selectedPhotons+process.egmGsfElectronIDSequence+process.egmPhotonIDSequence+process.photonIDValueMapProducer+process.electronIDValueMapProducer+process.photonMVAValueMapProducer+process.electronMVAValueMapProducer+process.JEC+process.ssbanalyzer)


