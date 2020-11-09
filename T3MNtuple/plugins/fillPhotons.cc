#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"

void T3MNtuple::fillPhotons(const edm::Event& iEvent,
                            const edm::EventSetup& iSetup,
                            const Handle<std::vector<reco::Photon> >& photons)
{

   for(unsigned int iPhoton = 0; iPhoton < photons->size(); iPhoton++){
      reco::PhotonRef photon(photons,iPhoton);
      std::vector<float>  iGammaP4;
      iGammaP4.push_back(photon->p4().E());
      iGammaP4.push_back(photon->p4().Px());
      iGammaP4.push_back(photon->p4().Py());
      iGammaP4.push_back(photon->p4().Pz());
      Gamma_P4.push_back(iGammaP4);

      Gamma_hasPixelSeed.push_back(photon->hasPixelSeed());
      Gamma_hasConversionTracks.push_back(photon->hasConversionTracks());
      Gamma_e1x5.push_back(photon->e1x5());
      Gamma_e2x5.push_back(photon->e2x5());
      Gamma_e3x3.push_back(photon->e3x3());
      Gamma_e5x5.push_back(photon->e5x5());

      Gamma_isPFPhoton.push_back(photon->isPFlowPhoton());

   }
}

void T3MNtuple::fillPhotons(const edm::Event& iEvent,
                            const edm::EventSetup& iSetup,
                            const Handle<std::vector<pat::Photon> >& photons)
{

   for(unsigned int iPhoton = 0; iPhoton < photons->size(); iPhoton++){
      pat::PhotonRef photon(photons,iPhoton);
      std::vector<float>  iGammaP4;
      iGammaP4.push_back(photon->p4().E());
      iGammaP4.push_back(photon->p4().Px());
      iGammaP4.push_back(photon->p4().Py());
      iGammaP4.push_back(photon->p4().Pz());
      Gamma_P4.push_back(iGammaP4);

      Gamma_hasPixelSeed.push_back(photon->hasPixelSeed());
      Gamma_hasConversionTracks.push_back(photon->hasConversionTracks());
      Gamma_e1x5.push_back(photon->e1x5());
      Gamma_e2x5.push_back(photon->e2x5());
      Gamma_e3x3.push_back(photon->e3x3());
      Gamma_e5x5.push_back(photon->e5x5());

      Gamma_isPFPhoton.push_back(photon->isPFlowPhoton());

   }
}
