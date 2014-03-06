import FWCore.ParameterSet.Config as cms

minBiasFiles_p1 = cms.untracked.vstring(
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/008E2E98-0A39-E311-833F-0025905938D4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/027029F2-FE38-E311-ACD6-003048678B34.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0298D3C4-9E3A-E311-BB62-003048B95B30.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/02D41CC7-9E3A-E311-979A-003048678FFA.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/02E67334-0639-E311-B292-003048679180.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/02ECE0B1-ED38-E311-A2F4-002618943916.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/02FF8ABA-9E3A-E311-B877-00304867924E.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/040D61A8-F438-E311-B316-003048FFCC18.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/04642097-FE38-E311-89CD-0025905938A4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/04FE79A1-0A39-E311-8670-002618943842.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/06C21BF5-FE38-E311-9653-00304867929E.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/06CA820D-F538-E311-9266-002590593920.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/06CD0E9A-0639-E311-96AF-0025905964B2.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/06DB478D-983A-E311-BD31-002618943843.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0802ADD9-EC38-E311-893F-0025905964C4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/085F4356-0439-E311-A044-003048FFD76E.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/08F0F097-0639-E311-858B-002618943882.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0AE5203A-0639-E311-9B69-003048FFCBA4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0C3AFEB2-483A-E311-A13D-002618943898.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0CCB409B-0639-E311-AF57-0025905964BC.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0E039F17-0B39-E311-8AB4-003048FFD7D4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0E078B9B-AA39-E311-85A2-002618FDA250.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0E1FF333-0639-E311-B58B-002618943947.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0E3A75F3-FE38-E311-A928-002618943924.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0E74B550-0439-E311-B5FA-00261894387E.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/100F8BBA-9E3A-E311-A3A6-00304867924E.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1046CEFB-0239-E311-9DE9-0025905964CC.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/105D57B9-9E3A-E311-AFEF-00304867BFBC.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/106DA150-0439-E311-95FC-002618943985.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/10A2F6B5-F138-E311-B229-003048FFCBFC.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/10A52FDD-963A-E311-A58C-002618943864.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/10B09A0A-F538-E311-A059-0026189438C1.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/12DDE416-0B39-E311-AE8B-002590596486.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/144AB833-0639-E311-8DA8-00261894396F.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/149ED172-F538-E311-B2CC-00261894380D.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/14AACA3A-0039-E311-8AFC-00261894382D.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/14D4FD93-FE38-E311-8F99-002618943959.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1674B648-0839-E311-AC1B-003048D3FC94.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/168AC4F1-FE38-E311-859F-002618943937.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/16A01D9F-0039-E311-9E13-002590596486.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/181DE50E-F538-E311-90F4-003048FFD756.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/189C427E-0339-E311-A6E0-003048678F9C.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/18E53376-F438-E311-8367-003048FFCBA8.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/18F7FD94-1739-E311-AF9E-00261894380D.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/18FEE8FD-0239-E311-9B5E-0025905964BC.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1A842285-0A39-E311-980F-0026189438A0.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1A9F2C3D-0039-E311-9261-002618943854.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1AD53F0E-F538-E311-BB0C-003048FF9AA6.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1C0CA584-0339-E311-84D6-003048FFCBA4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1C27DDF4-FE38-E311-B38B-00304867BED8.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1C34DE63-0439-E311-A230-003048678B04.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1C8501ED-0639-E311-A093-00261894383B.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1E7E3FA3-F438-E311-852C-002618943809.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2003DAAE-F438-E311-AB96-003048FFD7C2.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/22100BEE-353C-E311-9CF2-003048679076.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/22B45156-0439-E311-9984-003048FFD728.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/244564A0-0039-E311-ACDB-003048FFD7D4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/24BB9351-0439-E311-BDAC-00261894397D.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/24C56408-F538-E311-B143-0026189438C4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/24DEECA3-F438-E311-A474-0025905964C2.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/261B1AF1-FE38-E311-99DB-002618943950.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/280A107E-0339-E311-8855-00261894398C.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/28620798-FE38-E311-A6E8-003048FFD76E.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2A1D0C38-0639-E311-8FB9-0025905964BC.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2A6C8581-0339-E311-8B65-00248C55CC40.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2A7E7E75-F538-E311-9873-003048678B34.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2A9F1E97-FE38-E311-AB26-00304867D836.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2AE1AAED-0639-E311-A171-0026189438FA.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2C7CC036-0639-E311-A857-002590593876.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2E732C0B-F538-E311-BE3C-002618943935.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2EAB097E-0339-E311-BC48-00261894392D.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2EC316F4-FE38-E311-99F2-00261894386F.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/303F85F9-0239-E311-91A5-002618FDA279.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/30C6EC09-F538-E311-97FD-002618943885.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/30D6A375-F438-E311-BE68-003048FFCB84.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/30E24A33-0639-E311-B5DB-0026189438AB.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3202E4B3-ED38-E311-9AAB-0025905964C4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/320E6F36-0639-E311-B329-002590596498.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/325B4E53-0339-E311-897E-002590593878.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/32611EEB-0339-E311-83A1-003048FFD732.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/32B975E4-0A39-E311-824D-0025905964B4.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/342A82F5-353C-E311-9C5E-0026189438F2.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/363D47A6-F438-E311-A319-002590596468.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3679A598-FE38-E311-AECB-003048FFD71E.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/38150DA1-F538-E311-A8BD-003048FFCB96.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3817C10E-F538-E311-940E-003048FFCBB0.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3A9AD551-0439-E311-8981-002618943943.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3AEBF2A0-0039-E311-B0DD-003048FFD728.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3C4FDB97-0639-E311-A29E-002618943880.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3C8946A0-F438-E311-B768-00261894385D.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3C997151-0839-E311-A3B4-003048678B14.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3E2A71AE-9E3A-E311-ABE9-002618943958.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3EF102F5-453A-E311-BDFA-003048B835A2.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/400B7080-0539-E311-AE15-003048D3FC94.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/405A9B99-FE38-E311-8ACF-003048FFD756.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/40C180E0-963A-E311-9BF3-00259059642E.root",
"/store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/445AEBBC-9E3A-E311-BECC-0026189438FD.root"
)
