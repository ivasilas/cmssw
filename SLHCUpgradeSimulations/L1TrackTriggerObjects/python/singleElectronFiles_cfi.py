import FWCore.ParameterSet.Config as cms

singleElectronFiles = cms.untracked.vstring(
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/00D6C34E-0339-E311-836A-002618943880.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/06321333-FF38-E311-BE39-003048678ED2.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/06740333-FF38-E311-8C20-002618943945.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0852AA33-FF38-E311-9BBC-0026189438CE.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/08B3E078-FF38-E311-9FBD-00261894391C.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/08FE95D2-FF38-E311-BD72-003048FFD7BE.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0EDA0734-FF38-E311-8B96-00261894397D.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1238404F-0339-E311-AA48-003048D15E14.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/141AA1F5-FE38-E311-A5DE-003048FFCB96.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/18BC0A22-FF38-E311-9195-0026189438D2.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/18D9CE77-FF38-E311-A57C-0026189437FE.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1C4550B0-0739-E311-BBB1-003048678FD6.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1E131C24-973A-E311-8739-003048D15DB6.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/22F057F3-353C-E311-9204-0026189438E8.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2E9CD4F3-FE38-E311-937C-0026189438E1.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/30519134-FF38-E311-A46A-003048FFCC1E.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/32858933-FF38-E311-A273-002618943951.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/34B0BA34-FF38-E311-9CC1-003048FFD740.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3604A620-FF38-E311-BEF4-002618943862.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/36A92C50-0339-E311-81A4-002354EF3BDC.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/36BD9A50-0339-E311-8D2A-002618B27F8A.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/36D29CF7-FE38-E311-893B-003048FFD732.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3AA09A7C-FF38-E311-946A-003048FFD732.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3C1A3534-FF38-E311-9116-002618FDA216.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3CD3C5B2-9E3A-E311-BBC8-003048679164.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3E2B44B4-0739-E311-9733-003048678B94.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/42BE3778-FF38-E311-BD9C-00261894392F.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/44D53520-FF38-E311-841C-002618943900.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/48177BB4-0739-E311-8DCD-002590593902.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/4EA020F7-FE38-E311-912F-003048FF9AC6.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5028BE52-0339-E311-9963-002618943898.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/54111050-0339-E311-9CA7-003048678F62.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/54661A77-FF38-E311-9343-00261894386C.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/54E07033-FF38-E311-BB78-00261894394A.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5689A958-F638-E311-89EE-003048FFD71E.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/58366B34-FF38-E311-9FCB-003048FFCB96.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5A537579-FF38-E311-AC5B-003048679070.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5AC3CE78-FF38-E311-854C-00261894395C.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5CF5AC34-FF38-E311-AB71-003048FFCC1E.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5CF74477-FF38-E311-94E4-00261894382D.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5E420E34-FF38-E311-A859-0025905938D4.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5EF43933-FF38-E311-8307-0026189438E6.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5EFA53B0-0739-E311-B1A6-0030486791AA.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/60AC0C52-0339-E311-91AD-0026189438F6.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/6267FA32-FF38-E311-BDC5-00261894380B.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/64631277-FF38-E311-B82C-002618943905.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/66FDBF2A-973A-E311-BDFF-003048678F62.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/68702A05-973A-E311-B7B7-002590593878.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/6A1A3D21-FF38-E311-9EC8-00248C0BE018.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7068567C-FF38-E311-AFD1-003048FF9AA6.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/741F68F1-393A-E311-9733-00261894392F.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/742A1B33-FF38-E311-A02F-0026189437F5.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7668F432-FF38-E311-BDF7-002618FDA26D.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/786DF034-FF38-E311-8279-0026189438E7.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/78A2C633-FF38-E311-A65A-00304867906C.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/78F51578-FF38-E311-9E70-002618943826.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7C8E3255-0439-E311-AFDA-002590596486.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7CDBE478-FF38-E311-A919-002590596486.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/80381C08-453A-E311-9460-003048678ED2.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/82A46434-FF38-E311-8F92-003048FFCBA4.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/84B7CD33-FF38-E311-9171-002590596498.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/8E24CEB1-9E3A-E311-85FF-002618943932.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/8E2FC64F-0439-E311-AED6-002618943821.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/8E5C390F-FF38-E311-823E-003048FFD75C.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/90DB9C34-FF38-E311-A92F-003048FFD736.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/92B4FD7A-FF38-E311-9DF2-003048FFD71A.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/964FCEB2-0739-E311-9F59-0025905964C0.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9A12EC7B-FF38-E311-B1C9-003048FF86CA.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A24D2A0F-FF38-E311-9654-0025905964B2.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A6DDDEF3-FE38-E311-95C8-00304867C1BC.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/AA1924F6-FE38-E311-9387-003048D15D22.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/AE4C4E76-FF38-E311-8314-0026189437E8.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/AEAA6B7D-FF38-E311-95C2-00304867916E.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B0488834-FF38-E311-930D-003048FFD730.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B61E6A76-FF38-E311-9EDE-0026189438D7.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B6C3BE4F-0439-E311-B720-00304867C1BC.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B6DCF720-FF38-E311-BFD4-002618943833.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B6F2F921-FF38-E311-AF63-00261894387E.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BAF2DB55-0339-E311-B0AF-003048FFD71E.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BCB67F50-0339-E311-BCB9-002618943914.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BECC6533-FF38-E311-8912-00261894395B.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C6DF1BF4-FE38-E311-B210-002618943867.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C6FC9034-FF38-E311-98B1-003048FFD79C.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/CA64B134-FF38-E311-8070-003048FFCB84.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/CE5C1D50-0339-E311-A77F-003048678B86.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D061AE34-FF38-E311-9812-003048FFD796.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D8CE1950-0339-E311-AD95-003048678F9C.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/DE1BF350-0339-E311-8658-0026189438A9.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E0AADFB3-9E3A-E311-9ABD-002618943935.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E47DBA34-FF38-E311-8BB0-003048FFD796.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F4200834-FF38-E311-ADAE-0025905822B6.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F458B3F6-FE38-E311-A5E8-0026189438D5.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F6E656E6-0339-E311-A3DF-002590596484.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FE9BF034-FF38-E311-8584-0026189438E7.root",
"/store/mc/UpgFall13d/SingleElectronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FEDB1C0F-FF38-E311-A659-0025905938D4.root"
)

