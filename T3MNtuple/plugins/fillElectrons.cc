#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
#include "DsTau23Mu/T3MNtuple/interface/TrackParticle.h"
#include "DsTau23Mu/T3MNtuple/interface/ParticleBuilder.h"
#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"



void T3MNtuple::fillElectrons(const edm::Event& iEvent,
			      const edm::EventSetup& iSetup,
			      const Handle<TrackCollection>& trackCollection,
			      const Handle<VertexCollection>& pvs,
			      const Handle<BeamSpot>& beamSpotHandle,
			      const Handle<vector<pat::Electron> >& Electrons,
			      const Handle<vector<Vertex> >&  vertexs)
{
  unsigned int Electron_index = 0;
  for (vector<pat::Electron>::const_iterator iElectron = Electrons->begin(); iElectron != Electrons->end(); ++iElectron, Electron_index++) 
   {


     pat::ElectronRef RefElectron(Electrons, Electron_index);
     if(RefElectron->pt() > ElectronPtCut_)
       {
	 std::vector<float> iElectron_p4;
	 iElectron_p4.push_back(RefElectron->p4().E());
	 iElectron_p4.push_back(RefElectron->p4().Px());
	 iElectron_p4.push_back(RefElectron->p4().Py());
	 iElectron_p4.push_back(RefElectron->p4().Pz());
	 Electron_p4.push_back(iElectron_p4);
	 
	 Electron_Charge.push_back(RefElectron->charge());
	 // ----------------------   recomendations from 
	 // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_formats
	 
	 
	 Electron_puppiNeutralHadronIso.push_back(RefElectron->puppiNeutralHadronIso());
	 Electron_puppiChargedHadronIso.push_back(RefElectron->puppiChargedHadronIso());
	 Electron_puppiPhotonIso.push_back(RefElectron->puppiPhotonIso());
	 Electron_trackIso.push_back(RefElectron->trackIso());
	 Electron_isPF.push_back(RefElectron->isPF());
	 // combined isolation is taken from: 
	 // https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleEffAreaPFIsoCut.cc#L83-L94	 
	 //

	 // I ignore effective area correction, if needed look up here:
	 // https://indico.cern.ch/event/482673/contributions/2187022/attachments/1282446/1905912/talk_electron_ID_spring16.pdf
	 //

	 const reco::GsfElectron::PflowIsolationVariables& pfIso =
	   RefElectron->pfIsolationVariables();

	 //	 EffectiveAreas effectiveAreas;
	 //	 double absEta = std::abs(RefElectron->superCluster()->eta());
	 const float chad   = pfIso.sumChargedHadronPt;
	 const float nhad   = pfIso.sumNeutralHadronEt;
	 const float pho    = pfIso.sumPhotonEt;
	 //	 const float eA     = effectiveAreas.getEffectiveArea( absEta );
	 //	 const float rho    = _rhoHandle.isValid() ? (float)(*_rhoHandle) : 0; // std::max likes float arguments
 	 //	 const float iso    = chad + std::max(0.0f, nhad + pho - rho*eA);
	 const float iso     =  chad + std::max(0.0f, nhad + pho);
	 const float relEiso =  iso / (iso  + RefElectron->pt());
  
	 Electron_relativeIsolation.push_back(relEiso);
	 Electron_cutBasedElectronID_Fall17_94X_V2_veto.push_back(RefElectron->electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
	 Electron_cutBasedElectronID_Fall17_94X_V2_loose.push_back(RefElectron->electronID("cutBasedElectronID-Fall17-94X-V2-loose"));
	 Electron_cutBasedElectronID_Fall17_94X_V2_medium.push_back(RefElectron->electronID("cutBasedElectronID-Fall17-94X-V2-medium"));
	 Electron_cutBasedElectronID_Fall17_94X_V2_tight.push_back(RefElectron->electronID("cutBasedElectronID-Fall17-94X-V2-tight"));
       }
   }
}


