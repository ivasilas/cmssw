import FWCore.ParameterSet.Config as cms

process = cms.Process("Jet")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                            '/store/cmst3/user/eperez/L1TrackTrigger/612_SLHC6/muDST/TTbar/BE5D/zmatchingOff/m1_TTbar_BE5D.root'
                            )
                            )

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS261_V3::All', '')


# raw2digi to get the gct digis
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.p0 = cms.Path( process.RawToDigi )

# run L1Reco
process.load('Configuration.StandardSequences.L1Reco_cff')
process.L1Reco = cms.Path(process.l1extraParticles)


# run the L1TkJetProducer 
process.L1TkJets = cms.EDProducer("L1TkJetProducer",
                                  #L1CentralJetInputTag = cms.InputTag("l1extraParticles","Central"),
                                  L1CentralJetInputTag = cms.InputTag("L1JetsFromHIHLTJets",""),     
                                  L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
                                  TRK_ZMAX = cms.double(25.),     # max track z0 [cm]
                                  TRK_CHI2MAX = cms.double(100.), # max track chi2
                                  TRK_PTMIN = cms.double(2.0),    # minimum track pt [GeV]
                                  TRK_ETAMAX = cms.double(2.5),   # maximum track eta
                                  TRK_NSTUBMIN = cms.int32(4),    # minimum number of stubs on track
                                  TRK_NSTUBPSMIN = cms.int32(0),  # DOESN'T WORK, KEEP AT ZERO !!! minimum number of stubs in PS modules on track        
)
process.pJets = cms.Path(process.L1TkJets)


# run the L1TkHTMissProducer 
process.L1TkHTMiss = cms.EDProducer("L1TkHTMissProducer",
                                    L1TkJetInputTag = cms.InputTag("L1TkJets","Central"),
                                    JET_PTMIN = cms.double(20.0),          # minimum pt of jets for HT [GeV]
                                    JET_HLTETA = cms.bool(True),           # temporary hack, for HLT jets remove jets from bad eta regions
                                    DoVtxConstrain = cms.bool(True),       # turn on/off applying any vertex constraint
                                    DeltaZ = cms.double(1.0),              # require jets to have |z_jet - z_ref| < DeltaZ [cm]
                                    PrimaryVtxConstrain = cms.bool(False), # use primary vertex instead of leading jet as reference z position
                                    L1VertexInputTag = cms.InputTag("L1TrackPrimaryVertex")
)
process.pHTM = cms.Path(process.L1TkHTMiss)


process.Out = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "example.root" ),
    fastCloning = cms.untracked.bool( False ),
    outputCommands = cms.untracked.vstring( 'drop *')
)


process.Out.outputCommands.append('keep *_L1TkJets_*_*')
process.Out.outputCommands.append('keep *_L1TkHTMiss_*_*')

process.FEVToutput_step = cms.EndPath(process.Out)




