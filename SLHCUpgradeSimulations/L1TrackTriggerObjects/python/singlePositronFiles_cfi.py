import FWCore.ParameterSet.Config as cms

singlePositronFiles = cms.untracked.vstring(
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/00437C8C-0439-E311-8CF6-003048FFD760.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/009D288B-0439-E311-989D-0025905938D4.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/02C2173D-0039-E311-ADF0-002618943943.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/02CC0D84-0439-E311-9968-002354EF3BD0.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/040E54B3-9E3A-E311-AAB2-003048678FAE.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/060192BC-9E3A-E311-9DDE-00261894386F.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/06972787-0439-E311-9268-003048678FF6.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/06B03C86-0439-E311-9298-0026189437ED.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0E1F3189-0439-E311-86D2-002590593920.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/0EF0A282-0439-E311-9AC9-002618943933.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/103DD68C-0439-E311-9B00-003048FF9AC6.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/10DF35DD-973A-E311-83ED-003048FFCBA4.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/12E6BD86-0439-E311-B374-003048679188.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1AEE64C0-9E3A-E311-9217-003048678B44.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/1E7D0D85-0439-E311-987F-0026189438D6.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/246D0B87-0439-E311-931D-002618FDA204.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/24F23E88-0439-E311-8FBF-003048679162.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2E0AB5BD-9E3A-E311-A8BB-003048678ED4.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2E4663B0-0439-E311-A466-002354EF3BDC.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/2E52618D-0439-E311-958A-003048FFCC1E.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/301FE486-0439-E311-98B9-003048678B04.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/32F7AB8A-0439-E311-ABCE-002618943946.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3A1D11BB-9E3A-E311-891C-003048678BAE.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3C2D828D-0439-E311-9541-0025905938B4.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3C41588A-0439-E311-9B0D-002590593876.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/3CEE808C-0439-E311-9110-003048FFCC18.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/40329AC1-9E3A-E311-997E-003048678B94.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/405E0689-0439-E311-A6D5-002590593920.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/407A528A-0439-E311-B3AF-00259059391E.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/4698AB81-0339-E311-8C25-002618FDA248.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/482415B1-0439-E311-B64B-002618943831.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/50DAFB89-0439-E311-AB9B-00261894391B.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/52E0CB84-0439-E311-9F95-002618943939.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/58283441-0039-E311-9A79-003048FFCC18.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5868AF86-0439-E311-B254-003048FF9AC6.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/5A708EC0-9E3A-E311-88EE-002618943886.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/6CDADB37-403A-E311-BC19-00304867C1BA.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/746A72D2-3D3C-E311-9083-002618943986.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7647A485-0439-E311-8CE2-00261894394D.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7694FF87-0439-E311-BE4F-00261894390C.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/76BB21B2-9E3A-E311-A5C5-003048D15DB6.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7C66D084-0439-E311-A284-003048678BE6.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/7CF349DA-973A-E311-A1BC-002618943857.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/8A521983-0439-E311-8216-002618943979.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/90A2B9BE-9E3A-E311-ADE0-003048678EE2.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/94580186-0439-E311-9E5B-003048678C3A.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/94CE2087-0439-E311-B0F7-002618FDA265.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/96408883-0439-E311-9C1C-0026189438F3.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9834708C-0439-E311-9A54-003048FFCC2C.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/9CCEB88C-0439-E311-A27D-003048FFD770.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A07E36B0-0439-E311-A946-0026189438A9.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A251F085-0439-E311-94D2-003048679294.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A266CD84-0439-E311-B54C-00261894390E.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A2A71FB5-0439-E311-A532-0025905938A4.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A4780C87-0439-E311-8780-0026189438D3.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A4D7BEC0-9E3A-E311-B420-0030486792B8.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A61CA2B4-9E3A-E311-A808-002354EF3BE2.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A6E6568B-0439-E311-92A2-0026189438D9.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A854F78E-0439-E311-8034-002590593872.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A8661483-0439-E311-A056-0026189437F8.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A8C8C989-0439-E311-91D1-002590596498.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/A8D8A141-0039-E311-BFE1-003048FFD728.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/AA35D33D-0039-E311-BDD0-00261894393F.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/AA60A384-0439-E311-8207-0026189438C0.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B0316498-0839-E311-93DA-0025905964B6.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B40E598A-0439-E311-84EA-0026189437FC.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/B6EAD188-0439-E311-868F-002618943900.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BA3EC786-0439-E311-9DC1-002618943876.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BA5452B3-9E3A-E311-8D6D-0025905964C2.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BAA975B3-0439-E311-80A2-003048678BAC.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BE9FE5B2-9E3A-E311-9C53-003048678F62.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/BEDA8F87-0439-E311-B438-002618943807.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C2A0B73F-0039-E311-A60A-003048FFCB6A.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/C2E16D85-0439-E311-8B78-003048678B1A.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/CA111484-0439-E311-BAC7-00261894390E.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/CE300589-0439-E311-9E97-00259059642A.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D4B955CC-3D3C-E311-9DCD-0030486790B8.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D6BF4FB1-0439-E311-9068-002618943911.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D6BF5B87-0439-E311-A8CE-0026189438AA.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/D8E3BB88-0439-E311-9729-0030486790B8.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E06344B4-9E3A-E311-A763-00304867D838.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E2C37F85-0439-E311-A749-0025905964BA.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E899B07E-0339-E311-B996-0026189438D4.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/E8BECA84-0439-E311-943A-00261894397F.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/EA2BEF3B-0039-E311-B753-00261894394F.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/EA34D084-0439-E311-8B04-002618943809.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/EC263185-0439-E311-8C1E-002618943948.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/ECDFA28C-0439-E311-B744-003048FFCBA4.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/EE0C63B1-0439-E311-AFD4-003048678B38.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F0A3558A-0439-E311-B9D3-002590593876.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/F8233C8A-0439-E311-A43D-0026189437FC.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FAE22D87-0439-E311-B70D-002618943943.root",
"/store/mc/UpgFall13d/SinglePositronFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/FEED1685-0439-E311-BD96-00261894396D.root"
)
