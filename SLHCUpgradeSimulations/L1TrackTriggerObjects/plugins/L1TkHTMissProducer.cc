// -*- C++ -*-
//
//
// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Tue Nov 12 17:03:19 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticleFwd.h"


using namespace l1extra;

//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class L1TkHTMissProducer : public edm::EDProducer {
public:
  
  explicit L1TkHTMissProducer(const edm::ParameterSet&);
  ~L1TkHTMissProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  //virtual void beginRun(edm::Run&, edm::EventSetup const&);
  //virtual void endRun(edm::Run&, edm::EventSetup const&);
  //virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  //virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  
  // ---------- member data ---------------------------
  edm::InputTag L1TkJetInputTag;
  
  float JET_PTMIN;                // [GeV]
  float JET_ETAMAX;               // [rad]
  bool JET_HLTETA;                // temporary hack to remove bad jets when using HI HLT jets...
  bool DoVtxConstrain;            // require vertex constraint
  bool PrimaryVtxConstrain;       // use event primary vertex instead of leading jet (if DoVtxConstrain)
  edm::InputTag L1VertexInputTag; // used only when PrimaryVtxConstrain = True (if DoTvxConstrain)
  float DeltaZ;                   // for jets [cm] (if DoTvxConstrain)
  
};

//////////////
// constructor
L1TkHTMissProducer::L1TkHTMissProducer(const edm::ParameterSet& iConfig)
{
  L1TkJetInputTag = iConfig.getParameter<edm::InputTag>("L1TkJetInputTag");
  
  JET_PTMIN  = (float)iConfig.getParameter<double>("JET_PTMIN");
  JET_ETAMAX = (float)iConfig.getParameter<double>("JET_ETAMAX");
  JET_HLTETA = (float)iConfig.getParameter<bool>("JET_HLTETA");

  DoVtxConstrain      = iConfig.getParameter<bool>("DoVtxConstrain");
  PrimaryVtxConstrain = iConfig.getParameter<bool>("PrimaryVtxConstrain");
  L1VertexInputTag    = iConfig.getParameter<edm::InputTag>("L1VertexInputTag") ;
  DeltaZ              = (float)iConfig.getParameter<double>("DeltaZ");
  
  produces<L1TkHTMissParticleCollection>();
}

/////////////
// destructor
L1TkHTMissProducer::~L1TkHTMissProducer()
{
}



///////////
// producer
void L1TkHTMissProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
 
  // ----------------------------------------------------------------------------------------------
  // output container
  // ----------------------------------------------------------------------------------------------

  std::auto_ptr<L1TkHTMissParticleCollection> result(new L1TkHTMissParticleCollection);


  // ----------------------------------------------------------------------------------------------
  // retrieve input containers 
  // ----------------------------------------------------------------------------------------------

  // L1 primary vertex
  edm::Handle<L1TrackPrimaryVertexCollection> L1VertexHandle;
  iEvent.getByLabel(L1VertexInputTag,L1VertexHandle);
  std::vector<L1TrackPrimaryVertex>::const_iterator vtxIter;
  
  // L1 track-trigger jets
  edm::Handle<L1TkJetParticleCollection> L1TkJetsHandle;
  iEvent.getByLabel(L1TkJetInputTag,L1TkJetsHandle);
  std::vector<L1TkJetParticle>::const_iterator jetIter;
  

  if ( ! L1TkJetsHandle.isValid() ) {
    LogError("L1TkHTMissProducer")
      << "\nWarning: L1TkJetParticleCollection with " << L1TkJetInputTag
      << "\nrequested in configuration, but not found in the event. Exit"
      << std::endl;
    return;
  }
  

  // ----------------------------------------------------------------------------------------------
  // if PrimaryVtxConstrain, use the primary vertex instead of z position from leading jet
  // ----------------------------------------------------------------------------------------------

  float evt_zvtx = -999;
  bool found_vtx = false;

  edm::Ref< L1TrackPrimaryVertexCollection > L1VtxRef; 	// null reference

  if ( DoVtxConstrain && PrimaryVtxConstrain ) {
    if( !L1VertexHandle.isValid() ) {
      LogError("L1TkHTMissProducer")
	<< "\nWarning: L1TrackPrimaryVertexCollection with " << L1VertexInputTag
	<< "\nrequested in configuration, but not found in the event. Exit."
	<< std::endl;
      return ;
    }
    else {
      std::vector<L1TrackPrimaryVertex>::const_iterator vtxIter = L1VertexHandle->begin();
      // by convention, the first vertex in the collection is the one that should
      // be used by default
      evt_zvtx = vtxIter->getZvertex();
      found_vtx = true;
      int ivtx = 0;
      edm::Ref< L1TrackPrimaryVertexCollection > vtxRef( L1VertexHandle, ivtx );
      L1VtxRef = vtxRef;
    }
  } //endif PrimaryVtxConstrain


  // ----------------------------------------------------------------------------------------------
  // using z position of leading jet to define "event vertex"
  // ----------------------------------------------------------------------------------------------

  float zvtx_jetpt = -999;

  if ( DoVtxConstrain && !PrimaryVtxConstrain ) {

    for (jetIter = L1TkJetsHandle->begin(); jetIter != L1TkJetsHandle->end(); ++jetIter) {

      // only consider jets from the central BX
      int ibx = jetIter->bx();
      if (ibx != 0) continue;
      
      float jetVtx = jetIter->getJetVtx();
      float tmp_jet_pt  = jetIter->pt();
      float tmp_jet_eta = jetIter->eta();
      float tmp_jet_phi = jetIter->phi();
      
      if (tmp_jet_pt < JET_PTMIN) continue; 
      if (fabs(tmp_jet_eta) > JET_ETAMAX) continue;
      
      if (JET_HLTETA) {
	if ( (tmp_jet_eta > -2.0) && (tmp_jet_eta < -1.5) && (tmp_jet_phi > -1.7) && (tmp_jet_phi < -1.3) ) continue;
	if ( (tmp_jet_eta > 0.0) && (tmp_jet_eta < 1.4) && (tmp_jet_phi > -2.0) && (tmp_jet_phi < -1.6) ) continue;
      }
      
      // find vertex position of leading jet
      if (tmp_jet_pt > zvtx_jetpt) {
	evt_zvtx = jetVtx;
	zvtx_jetpt = tmp_jet_pt;
	found_vtx = true; 
      }

    }//end loop over jets

  }//endif z position from leading jet


  float sumPx = 0;
  float sumPy = 0;
  float etTot = 0; //HT
  

  if (!found_vtx) std::cout << "WARNING: didn't find any z vertex for this event!" << std::endl;  



  // ----------------------------------------------------------------------------------------------
  // loop over jets
  // ----------------------------------------------------------------------------------------------

  for (jetIter = L1TkJetsHandle->begin(); jetIter != L1TkJetsHandle->end(); ++jetIter) {

    // only consider jets from the central BX
    int ibx = jetIter->bx();
    if (ibx != 0) continue;
    
    float px = jetIter->px();
    float py = jetIter->py();
    float et = jetIter->et();
    float jetVtx = jetIter->getJetVtx();
    float tmp_jet_pt  = jetIter->pt();
    float tmp_jet_eta = jetIter->eta();
    float tmp_jet_phi = jetIter->phi();
    
    if (tmp_jet_pt < JET_PTMIN) continue; 
    if (fabs(tmp_jet_eta) > JET_ETAMAX) continue;

    if (JET_HLTETA) {
      if ( (tmp_jet_eta > -2.0) && (tmp_jet_eta < -1.5) && (tmp_jet_phi > -1.7) && (tmp_jet_phi < -1.3) ) continue;
      if ( (tmp_jet_eta > 0.0) && (tmp_jet_eta < 1.4) && (tmp_jet_phi > -2.0) && (tmp_jet_phi < -1.6) ) continue;
    }

    // vertex consistency requirement
    bool VtxRequirement = fabs(jetVtx - evt_zvtx) < DeltaZ;
    
    if (!DoVtxConstrain || VtxRequirement) {
      sumPx += px;
      sumPy += py;
      etTot += et;
    }
    
  }//end loop over jets


  // ----------------------------------------------------------------------------------------------
  // define missing HT 
  // ----------------------------------------------------------------------------------------------

  float et = sqrt(sumPx*sumPx + sumPy*sumPy);
  math::XYZTLorentzVector missingEt(-sumPx, -sumPy, 0, et); 

  edm::RefProd<L1TkJetParticleCollection> jetCollRef(L1TkJetsHandle);
  L1TkHTMissParticle tkHTM(missingEt, etTot, jetCollRef, L1VtxRef);

  if (DoVtxConstrain && !PrimaryVtxConstrain) {
    tkHTM.setVtx(evt_zvtx);
  }

  result->push_back(tkHTM);
  iEvent.put(result);

}


// ------------ method called once each job just before starting event loop  ------------
void 
L1TkHTMissProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TkHTMissProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TkHTMissProducer::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
L1TkHTMissProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TkHTMissProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TkHTMissProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TkHTMissProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TkHTMissProducer);
