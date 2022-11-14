#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
#include "DsTau23Mu/T3MNtuple/interface/TrackParticle.h"
#include "DsTau23Mu/T3MNtuple/interface/ParticleBuilder.h"

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
	 
	 
	 Electron_cutBasedElectronID_Fall17_94X_V2_veto.push_back(RefElectron->electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
	 Electron_cutBasedElectronID_Fall17_94X_V2_loose.push_back(RefElectron->electronID("cutBasedElectronID-Fall17-94X-V2-loose"));
	 Electron_cutBasedElectronID_Fall17_94X_V2_medium.push_back(RefElectron->electronID("cutBasedElectronID-Fall17-94X-V2-medium"));
	 Electron_cutBasedElectronID_Fall17_94X_V2_tight.push_back(RefElectron->electronID("cutBasedElectronID-Fall17-94X-V2-tight"));
       }
   }
}


