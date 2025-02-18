'''
Based on the "reHLT" menu:
/frozen/2016/25ns10e33/v2.1/HLT/V3 (CMSSW_8_0_11, HLTrigger/Configuration/python/HLT_25ns10e33_v2_cff.py)
link: https://cmsweb-testbed.cern.ch/confdb/#config=/frozen/2016/25ns10e33/v2.1/HLT/V3
'''

triggerTable = {
    "ZnnHbb" : [
        "HLT_PFMET90_PFMHT90_IDTight_v*",
        "HLT_PFMET170_NoiseCleaned_v*",
    ],
    "ZnnHbbAll" : [
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_v*",
        "HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_v*",
        "HLT_PFMET90_PFMHT90_IDTight_v*",
        "HLT_PFMET100_PFMHT100_IDTight_v*",
        "HLT_PFMET110_PFMHT110_IDTight_v*",
        "HLT_PFMET120_PFMHT120_IDTight_v*",
        "HLT_PFMET170_NoiseCleaned_v*",
        "HLT_PFMET170_HBHE_BeamHaloCleaned_v*",
        "HLT_PFMET170_HBHECleaned_v*",
        "HLT_DiCentralPFJet55_PFMET110_v*",
        "HLT_PFHT350_PFMET100_v*",
    ],
    "ZeeHbbAll" : [
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Ele32_eta2p1_WPTight_Gsf_v*",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v*",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
        "HLT_Ele27_eta2p1_WPTight_Gsf_v*",
        "HLT_Ele22_eta2p1_WPLoose_Gsf_v*",
        "HLT_Ele23_WPLoose_Gsf_v*",
        "HLT_Ele27_WPLoose_Gsf_v*",
        "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
        "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*",
        "HLT_Ele30WP60_Ele8_Mass55_v*"
    ],
    "ZeeHbbHighLumi" : [
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*"
    ],
    "ZeeHbbLowLumi" : [
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*"
    ],


    "WenHbbAll" : [
        "HLT_Ele32_eta2p1_WPTight_Gsf_v*",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v*",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
        "HLT_Ele27_eta2p1_WPTight_Gsf_v*",
        "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
        "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*",
        "HLT_Ele27_WPLoose_Gsf_WHbbBoost_v*",
        "HLT_Ele22_eta2p1_WPLoose_Gsf_v*",
        "HLT_Ele23_WPLoose_Gsf_v*",
        "HLT_Ele27_WPLoose_Gsf_v*",
        "HLT_Ele27_WPTight_Gsf_v*",
        "HLT_Ele23_WPLoose_Gsf_WHbbBoost_v*",
        "HLT_Ele25_WPTight_Gsf_v*",
        "HLT_Ele25_eta2p1_WPLoose_Gsf_v*",
        "HLT_Ele25_eta2p1_WPTight_Gsf_v*" 
    ],
    "WenHbbHighLumi" : [
        "HLT_Ele27_WPTight_Gsf_v*",
        "HLT_Ele25_eta2p1_WPTight_Gsf_v*" ,
        "HLT_Ele32_eta2p1_WPTight_Gsf_v*",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v*",
        "HLT_Ele27_WPLoose_Gsf_WHbbBoost_v*",
    ],
    "WenHbbLowLumi" : [
        "HLT_Ele27_eta2p1_WPTight_Gsf_v*",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v*",
    ],

    "ZmmHbbAll" : [
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TkMu8_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
        "HLT_DoubleIsoMu17_eta2p1_v*",
        "HLT_IsoMu24_v*",
        "HLT_Mu24_eta2p1_v*",
        "HLT_TkMu24_eta2p1_v*",
        "HLT_IsoMu27_v*",
        "HLT_IsoTkMu27_v*",
        "HLT_IsoTkMu24_v*",
        "HLT_TkMu27_v*",
        "HLT_Mu27_v*",
        "HLT_IsoMu20_v*",
        "HLT_Mu20_v*",
        "HLT_TkMu20_v*",
        "HLT_IsoMu20_v*",
        "HLT_IsoTkMu20_v*",
        "HLT_Mu20_v*",
        "HLT_TkMu20_v*",
        "HLT_Mu40_eta2p1_PFJet200_PFJet50_v*",
        "HLT_IsoMu18_v*",
        "HLT_IsoTkMu18_v*",
    ],
    "ZmmHbbHighLumi" : [
        "HLT_IsoMu24_v*",
        "HLT_IsoMu24_eta2p1_v*",
        "HLT_IsoMu27_v*",
        "HLT_IsoTkMu27_v*",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
    ],
    "ZmmHbbLowLumi" : [
        "HLT_IsoMu20_v*",
        "HLT_IsoTkMu20_v*",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
    ],

    "WmnHbbAll" : [
        "HLT_Mu24_eta2p1_v*",
        "HLT_TkMu24_eta2p1_v*",
        "HLT_IsoMu27_v*",
        "HLT_IsoTkMu27_v*",
        "HLT_IsoTkMu24_v*",
        "HLT_TkMu27_v*",
        "HLT_Mu27_v*",
        "HLT_IsoMu20_v*",
        "HLT_Mu20_v*",
        "HLT_TkMu20_v*",
        "HLT_IsoMu20_v*",
        "HLT_IsoTkMu20_v*",
        "HLT_Mu20_v*",
        "HLT_TkMu20_v*",
        "HLT_IsoMu18_v*",
        "HLT_IsoTkMu18_v*",
        "HLT_Mu40_eta2p1_PFJet200_PFJet50_v*",
        "HLT_IsoMu16_eta2p1_MET30_v*",
        "HLT_Mu16_eta2p1_MET30_v*",
        "HLT_PFMET120_Mu5_v*",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",
        "HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight_v*",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*",
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v*",
        "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v*",
        "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v*",
    ],
    "WmnHbbHighLumi" : [
        "HLT_IsoMu24_v*",
        "HLT_IsoMu24_eta2p1_v*",
        "HLT_IsoMu27_v*",
        "HLT_IsoTkMu27_v*",
        "HLT_IsoTkMu24_v*",
        "HLT_IsoMu16_eta2p1_MET30_v*"
    ],
    "WmnHbbLowLumi" : [
        "HLT_IsoMu20_v*",
        "HLT_IsoTkMu20_v*",
        "HLT_IsoMu16_eta2p1_MET30_v*",
    ],

    "WtaunHbbAll" : [
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v*",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v*",
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_v*",
    ],
    "WtaunHbbHighLumi" : [
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v*",
    ],
    "WtaunHbbLowLumi" : [
        "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80_v*",
    ],

    "VBFHbbAll" : [
        "HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_v*",
        "HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200_v*",
        "HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_v*",
        "HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460_v*",
        "HLT_QuadPFJet_VBF_v*",
        "HLT_L1_TripleJet_VBF_v*",
        "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*",
        "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v*",
    ],
    "VBFHbbHighLumi" : [
        "HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_v*",
        "HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_v*",
    ],
    "VBFHbbLowLumi" : [
        "HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200_v*",
        "HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460_v*",
    ],

    "HH4bAll" : [
        "HLT_QuadJet45_TripleBTagCSV_p087_v*",
        "HLT_QuadJet45_DoubleBTagCSV_p087_v*",
        "HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v*",
        "HLT_DoubleJet90_Double30_DoubleBTagCSV_p087_v*",
        "HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v*",
        "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v*",
        "HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v*",
        "HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v*",
        "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v*",
        "HLT_AK8PFJet360_TrimMass30_v*",
        "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*",
        "HLT_PFHT800_v*",
        "HLT_PFHT900_v*"
    ],
    "HH4bHighLumi" : [
        "HLT_QuadJet45_TripleBTagCSV_p087_v*",
        "HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v*",
    ],
    "HH4bLowLumi" : [
        "HLT_QuadJet45_TripleBTagCSV_p087_v*",
        "HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v*",
    ],
    "ttH_SL_el" : [
        "HLT_Ele27_eta2p1_WPTight_Gsf_v*",
    ],
    "ttH_SL_mu" : [
        "HLT_IsoMu22_v*",
        "HLT_IsoTkMu22_v*",
    ],
    "ttH_DL_mumu" : [
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
    ],
    "ttH_DL_elmu" : [
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
    ],
    "ttH_DL_elel" : [
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
    ],
    "ttH_FH" : [
        "HLT_PFHT450_SixJet40_BTagCSV_p056_v*",
        "HLT_PFHT400_SixJet30_DoubleBTagCSV_p056_v*",
    ],
    "ttH_FH_prescaled" : [
        "HLT_PFHT450_SixJet40_v*",
        "HLT_PFHT400_SixJet30_v*",
    ],

    "ttH_htt" : [
        "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v*",
        "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v*",
        "HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v*",
        "HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v*",
        "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v*",
        "HLT_VLooseIsoPFTau120_eta2p1_v*",
        "HLT_VLooseIsoPFTau140_eta2p1_v*",
    ],
    "ttH_htt_lepFakeBgr" : [
        "HLT_Mu3_PFJet40_v*",
        "HLT_Mu8_v*",
        "HLT_Mu17_v*",
        "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*",
        "HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v*",
        "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*",
        "HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v*"
    ],

    "hadronic" : [
        "HLT_PFHT750_4JetPt50_v*",
        "HLT_PFHT800_v*",
        "HLT_PFHT900_v*",
        "HLT_PFJet40_v*",
        "HLT_PFJet60_v*",
        "HLT_PFJet80_v*",
        "HLT_PFJet140_v*",
        "HLT_PFJet200_v*",
        "HLT_PFJet260_v*",
        "HLT_PFJet320_v*",
        "HLT_PFJet400_v*",
        "HLT_PFJet450_v*",

        "HLT_DiPFJetAve40_v*",
        "HLT_DiPFJetAve60_v*",
        "HLT_DiPFJetAve80_v*",
        "HLT_DiPFJetAve140_v*",
        "HLT_DiPFJetAve200_v*",
        "HLT_DiPFJetAve260_v*",
        "HLT_DiPFJetAve320_v*",
        "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_v*",
        "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v*"

    ],
}
