import FWCore.ParameterSet.Config as cms

genSimAnalyzer = cms.EDAnalyzer('genSimAnalyzer',
                                fastTimeBarrelHits  = cms.InputTag("", "", ""), 
                                fastTimeEndcapHits  = cms.InputTag("", "", ""), 
                                genParticles        = cms.InputTag("", "", ""),
                                genParticles_t      = cms.InputTag("", "", ""),
                                simTracks           = cms.InputTag("", "", ""),
                                pdgIdTest           = cms.int32(13)
                                )
