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
    std::cout<<"  Electron index "<< Electron_index << " pT   "<< RefElectron->p4().E() << std::endl;

    }








}


