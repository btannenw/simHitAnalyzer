import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('L1',eras.Phase2_timing_layer_new)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# for MTD geometry
process.load("Geometry.CMSCommonData.cmsRecoIdealGeometryXML_cfi")
#process.load("Geometry.MTDNumberingBuilder.MTDModuleNumbering_cfi")
process.load("Geometry.MTDNumberingBuilder.mtdNumberingGeometry_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdGeometry_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        #'root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValSingleMuFlatPt_0p7to10_pythia8/GEN-SIM/103X_upgrade2023_realistic_v2_2023D35noPU-v1/10000/F4E8C5B3-3385-474F-8C40-8657D0E22F23.root' # muon
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_4_0_pre2/RelValSingleKaonFlatPt_0p7to10_pythia8_cfi/GEN-SIM/103X_upgrade2023_realistic_v2_2023D35noPU-v1/10000/FB6D9EFA-7D67-2444-829F-50E95B2046C7.root' # kaon
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_4_0_pre2/RelValSingleProtonFlatPt_0p7to10_pythia8_cfi/GEN-SIM/103X_upgrade2023_realistic_v2_2023D35noPU-v1/10000/9399E811-002E-8D48-8A45-545FA6D6954C.root' #proton
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_4_0_pre2/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/GEN-SIM/103X_upgrade2023_realistic_v2_2023D35noPU-v1/10000/D80EC4E5-0172-B648-A0AD-E9FCCEA66639.root' # pion
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_4_0_pre2/RelValSingleMuFlatPt_0p7to10_pythia8/GEN-SIM-DIGI-RAW/103X_upgrade2023_realistic_v2_2023D35noPU-v1/10000/4A5CE6E4-D298-4E4C-8601-36C87D732A99.root' # GEN-SIM-DIGI-RAW, muon
        ),
                            inputCommands = cms.untracked.vstring("keep *", 
        "drop l1tHGCalTowerMapBXVector_hgcalTriggerPrimitiveDigiProducer_towerMap_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT",
        "drop l1tEMTFHit2016s_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016s_simEmtfDigis__HLT")
     #skipEvents = cms.untracked.uint32(80)
                            )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives)


##the analyzer
process.genSimAnalyzer = cms.EDAnalyzer('genSimAnalyzer',
                                       fastTimeBarrelHits  = cms.InputTag("g4SimHits", "FastTimerHitsBarrel", "SIM"), 
                                       fastTimeEndcapHits  = cms.InputTag("g4SimHits", "FastTimerHitsEndcap", "SIM"), 
                                       genParticles        = cms.InputTag("genParticles", "", "SIM"),
                                       genParticles_t      = cms.InputTag("genParticles", "t0", "SIM"),
                                       simTracks           = cms.InputTag("g4SimHits", "", "SIM"),
                                       pdgIdTest           = cms.int32(13)
                                       )

process.analyzer = cms.Path(process.genSimAnalyzer)


process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("analyzer_singleProton_modType.root"), 
                                   closeFileFast = cms.untracked.bool(True)
)

# Schedule definition
#process.schedule = cms.Schedule(process.EcalEBtp_step,process.L1TrackTrigger_step,process.L1simulation_step,process.timingtracks,process.l1pf,process.L1PFTaus,process.analyzer,process.endjob_step) 
process.schedule = cms.Schedule(process.analyzer) 

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion



