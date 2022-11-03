#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
#include "DsTau23Mu/T3MNtuple/interface/LorentzVectorParticle.h"
#include "TLorentzVector.h"

int T3MNtuple::fillThreeMuons(const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<vector<reco::Muon>>& muons,
      const Handle<TrackCollection>& tracks,
      const Handle<BeamSpot>& beamSpotHandle,
      const Handle<trigger::TriggerEvent>& triggerSummary)
{

   std::vector<std::vector<unsigned int> > PreselectedThreeMuonsCollection = findThreeMuonsCandidates(iEvent, iSetup, beamSpotHandle, muons, trackCollection);
   if(PreselectedThreeMuonsCollection.size()==0){
      return 0;            //No three muons candidate found! Skip the event
   }

   for ( auto &iThreeMuon :  PreselectedThreeMuonsCollection ) {
      vector<TransientTrack> t_trks;   
      TransientVertex transVtx;
      ESHandle<TransientTrackBuilder> theB;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
      reco::MuonRef Muon1(muons, iThreeMuon.at(0));
      reco::MuonRef Muon2(muons, iThreeMuon.at(1));
      reco::MuonRef Muon3(muons, iThreeMuon.at(2));

      TrackRef track1 = Muon1->innerTrack();
      TrackRef track2 = Muon2->innerTrack();
      TrackRef track3 = Muon3->innerTrack();

      t_trks.push_back(theB->build(track1));
      t_trks.push_back(theB->build(track2));
      t_trks.push_back(theB->build(track3));

      KalmanVertexFitter kvf(true);
      bool FitOk(true);
      try {
         transVtx = kvf.vertex(t_trks); //KalmanVertexFitter
      } catch (...) {
         FitOk = false;
      }
      if (!transVtx.hasRefittedTracks())
         FitOk = false;
      if (transVtx.refittedTracks().size() != t_trks.size())
         FitOk = false;

      if(transVtx.isValid()){
         if(transVtx.totalChiSquared() < 100.)   //remove for sync
         { // very loose/ ndf =3

            int ntp = signalTau_lvp.size();
            signalTau_lvp.push_back(std::vector<double>());
            signalTau_cov.push_back(std::vector<double>());
            if(FitOk){
               signalTau_isLVP.push_back(1);
               LorentzVectorParticle  signalTau;
               GlobalPoint sv(transVtx.position().x(), transVtx.position().y(), transVtx.position().z());
               KinematicParticleFactoryFromTransientTrack kinFactory;
               float muMassSigma(sqrt(pow(10., -12.))), piChi(0.0), piNdf(0.0);
               std::vector<RefCountedKinematicParticle> muons;
               for (unsigned int i = 0; i <t_trks.size(); i++)
                  muons.push_back(kinFactory.particle(t_trks.at(i), PDGInfo::mu_mass(), piChi, piNdf, sv, muMassSigma));
               KinematicParticleVertexFitter kpvFitter;
               RefCountedKinematicTree jpTree = kpvFitter.fit(muons);
               if(jpTree->isValid()){
                  jpTree->movePointerToTheTop();
                  const KinematicParameters parameters = jpTree->currentParticle()->currentState().kinematicParameters();
                  AlgebraicSymMatrix77 cov = jpTree->currentParticle()->currentState().kinematicParametersError().matrix();

                  double c(0);
                  for (unsigned int i = 0; i < t_trks.size(); i++) {
                     c += t_trks.at(i).charge();
                  }

                  TMatrixT<double> tau_par(LorentzVectorParticle::NLorentzandVertexPar, 1);
                  TMatrixTSym<double> tau_cov(LorentzVectorParticle::NLorentzandVertexPar);
                  for (int i = 0; i < LorentzVectorParticle::NLorentzandVertexPar; i++) {
                     tau_par(i, 0) = parameters(i);
                     for (int j = 0; j < LorentzVectorParticle::NLorentzandVertexPar; j++) {
                        tau_cov(i, j) = cov(i, j);
                     }
                  }
                  signalTau = LorentzVectorParticle(tau_par, tau_cov, abs(PDGInfo::tau_minus) * c, c, theB->field()->inInverseGeV(sv).z());
                  signalTau_charge.push_back(signalTau.Charge());
                  signalTau_pdgid.push_back(signalTau.PDGID());
                  signalTau_B.push_back(signalTau.BField());
                  signalTau_M.push_back(signalTau.Mass());

                  for (int i = 0; i < signalTau.NParameters(); i++) {
                     signalTau_lvp.at(ntp).push_back(signalTau.Parameter(i));
                     for (int j = i; j < signalTau.NParameters(); j++) {
                        //	      signalTau_cov.at(ntp).push_back(signalTau.Covariance(i, j));// comment out to keep size low
                     }
                  }
               }
            }else{ signalTau_isLVP.push_back(-1);}

            ThreeMuons_idx.push_back(iThreeMuon);
            ThreeMuons_SV_Chi2.push_back(transVtx.totalChiSquared());
            ThreeMuons_SV_NDF.push_back(transVtx.degreesOfFreedom());

            std::vector<float> iTrigMatchdR(3, 999.0);
            vector<TLorentzVector> MuonTriggMatch;

            for ( auto &iMuon :  iThreeMuon ) {
               reco::MuonRef tmpRef(muons, iMuon);
               double tmpRef_p = tmpRef->p();
               TLorentzVector tmpP4(tmpRef->px(), tmpRef->py(), tmpRef->pz(),
                                    sqrt(pow(tmpRef_p,2)+pow(PDGInfo::mu_mass(),2)));
               MuonTriggMatch.push_back(tmpP4);
               TLorentzVector mut;
            }
            TriggerMatch(triggerSummary,  MuonTriggMatch, TriggerMuonMatchingdr_, iTrigMatchdR);
            ThreeMuons_TriggerMatch_dR.push_back(iTrigMatchdR);
         }
      }
   }
   return ThreeMuons_idx.size();
}

int T3MNtuple::fillThreeMuons(const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<vector<pat::Muon>>& muons,
      const Handle<TrackCollection>& tracks,
      const Handle<BeamSpot>& beamSpotHandle,
      const Handle<vector<pat::TriggerObjectStandAlone>>& triggerObjects,
      const TriggerNames& triggerNames)
{

   std::vector<std::vector<unsigned int> > PreselectedThreeMuonsCollection = findThreeMuonsCandidates(iEvent, iSetup, beamSpotHandle, muons, trackCollection);
   if (DEBUG) cout<<"preselected 3mu size = "<<PreselectedThreeMuonsCollection.size()<<endl;
   if(PreselectedThreeMuonsCollection.size()==0){
      return 0;            //No three muons candidate found! Skip the event
   }

   for ( auto &iThreeMuon :  PreselectedThreeMuonsCollection ) {
      vector<TransientTrack> t_trks;   
      TransientVertex transVtx;
      ESHandle<TransientTrackBuilder> theB;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
      pat::MuonRef Muon1(muons, iThreeMuon.at(0));
      pat::MuonRef Muon2(muons, iThreeMuon.at(1));
      pat::MuonRef Muon3(muons, iThreeMuon.at(2));

      TrackRef track1 = Muon1->innerTrack();
      TrackRef track2 = Muon2->innerTrack();
      TrackRef track3 = Muon3->innerTrack();

      t_trks.push_back(theB->build(track1));
      t_trks.push_back(theB->build(track2));
      t_trks.push_back(theB->build(track3));

      KalmanVertexFitter kvf(true);
      bool FitOk(true);
      try {
         transVtx = kvf.vertex(t_trks); //KalmanVertexFitter
      } catch (...) {
         FitOk = false;
      }
      if (!transVtx.hasRefittedTracks())
         FitOk = false;
      if (transVtx.refittedTracks().size() != t_trks.size())
         FitOk = false;

      if(transVtx.isValid()){
         if(transVtx.totalChiSquared() < 100.)   //remove for sync
         { // very loose/ ndf =3

            int ntp = signalTau_lvp.size();
            signalTau_lvp.push_back(std::vector<double>());
            signalTau_cov.push_back(std::vector<double>());
            if(FitOk){
               signalTau_isLVP.push_back(1);
               LorentzVectorParticle  signalTau;
               GlobalPoint sv(transVtx.position().x(), transVtx.position().y(), transVtx.position().z());
               KinematicParticleFactoryFromTransientTrack kinFactory;
               float muMassSigma(sqrt(pow(10., -12.))), piChi(0.0), piNdf(0.0);
               std::vector<RefCountedKinematicParticle> muons;
               for (unsigned int i = 0; i <t_trks.size(); i++)
                  muons.push_back(kinFactory.particle(t_trks.at(i), PDGInfo::mu_mass(), piChi, piNdf, sv, muMassSigma));
               KinematicParticleVertexFitter kpvFitter;
               RefCountedKinematicTree jpTree = kpvFitter.fit(muons);
               if(jpTree->isValid()){
                  jpTree->movePointerToTheTop();
                  const KinematicParameters parameters = jpTree->currentParticle()->currentState().kinematicParameters();
                  AlgebraicSymMatrix77 cov = jpTree->currentParticle()->currentState().kinematicParametersError().matrix();

                  double c(0);
                  for (unsigned int i = 0; i < t_trks.size(); i++) {
                     c += t_trks.at(i).charge();
                  }

                  TMatrixT<double> tau_par(LorentzVectorParticle::NLorentzandVertexPar, 1);
                  TMatrixTSym<double> tau_cov(LorentzVectorParticle::NLorentzandVertexPar);
                  for (int i = 0; i < LorentzVectorParticle::NLorentzandVertexPar; i++) {
                     tau_par(i, 0) = parameters(i);
                     for (int j = 0; j < LorentzVectorParticle::NLorentzandVertexPar; j++) {
                        tau_cov(i, j) = cov(i, j);
                     }
                  }
                  signalTau = LorentzVectorParticle(tau_par, tau_cov, abs(PDGInfo::tau_minus) * c, c, theB->field()->inInverseGeV(sv).z());
                  signalTau_charge.push_back(signalTau.Charge());
                  signalTau_pdgid.push_back(signalTau.PDGID());
                  signalTau_B.push_back(signalTau.BField());
                  signalTau_M.push_back(signalTau.Mass());

                  for (int i = 0; i < signalTau.NParameters(); i++) {
                     signalTau_lvp.at(ntp).push_back(signalTau.Parameter(i));
                     for (int j = i; j < signalTau.NParameters(); j++) {
		       signalTau_cov.at(ntp).push_back(signalTau.Covariance(i, j));// comment out to keep size low
                     }
                  }
               }
            }else{ signalTau_isLVP.push_back(-1);}

            ThreeMuons_idx.push_back(iThreeMuon);
            ThreeMuons_SV_Chi2.push_back(transVtx.totalChiSquared());
            ThreeMuons_SV_NDF.push_back(transVtx.degreesOfFreedom());

            std::vector<float> iTrigMatchdR(3, 999.0);
            vector<TLorentzVector> MuonTriggMatch;

            for ( auto &iMuon :  iThreeMuon ) {
               pat::MuonRef tmpRef(muons, iMuon);
               double tmpRef_p = tmpRef->p();
               TLorentzVector tmpP4(tmpRef->px(), tmpRef->py(), tmpRef->pz(),
                                    sqrt(pow(tmpRef_p,2)+pow(PDGInfo::mu_mass(),2)));
               MuonTriggMatch.push_back(tmpP4);
               TLorentzVector mut;
            }
            TriggerMatch(iEvent, triggerObjects, triggerNames , MuonTriggMatch, TriggerMuonMatchingdr_, iTrigMatchdR);
            ThreeMuons_TriggerMatch_dR.push_back(iTrigMatchdR);
         }
      }
   }
   return ThreeMuons_idx.size();
}
