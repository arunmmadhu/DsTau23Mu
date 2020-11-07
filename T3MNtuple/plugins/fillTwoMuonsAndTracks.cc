#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
#include "TLorentzVector.h"

int T3MNtuple::fillTwoMuonsAndTracks(const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<vector<reco::Muon> >& muons,
      const Handle<TrackCollection>& tracks,
      const Handle<BeamSpot>& beamSpotHandle,
      const Handle<trigger::TriggerEvent>& triggerSummary)
{

   std::vector<std::vector<unsigned int> > PreselectedTwoMuonsTrackCollection  = findTwoMuonsAndTrackCandidates(iEvent, iSetup, beamSpotHandle, muons, trackCollection);
   if(DEBUG)  std::cout<<" PreselectedTwoMuonsTrackCollection  "<< PreselectedTwoMuonsTrackCollection.size()<< std::endl;
   if(PreselectedTwoMuonsTrackCollection.size()==0){
      return 0;          //No two muons + track candidate found!
   }



   for ( auto &iTwoMuTr :  PreselectedTwoMuonsTrackCollection ) {
      vector<TransientTrack> t_trks;
      TransientVertex transVtx;
      ESHandle<TransientTrackBuilder> theB;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
      reco::MuonRef Muon1(muons, iTwoMuTr.at(0));
      reco::MuonRef Muon2(muons, iTwoMuTr.at(1));

      TLorentzVector mv1,mv2;

      mv1.SetPtEtaPhiM(Muon1->pt(), Muon1->eta(), Muon1->phi(), 0.106);
      mv2.SetPtEtaPhiM(Muon2->pt(), Muon2->eta(), Muon2->phi(), 0.106);

      TrackRef track1 = Muon1->innerTrack();
      TrackRef track2 = Muon2->innerTrack();
      TrackRef track3 = TrackRef(trackCollection, iTwoMuTr.at(2));

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



      if(FitOk){
         if(transVtx.totalChiSquared() < 100){
            if( (mv1 + mv2).M() < phimassmin_)// || (mv1 + mv2).M() > phimassmax_)   // renmove that
            {
               TwoMuonsTrack_idx.push_back(iTwoMuTr);
               TwoMuonsTrack_SV_Chi2.push_back(transVtx.totalChiSquared());
               TwoMuonsTrack_SV_NDF.push_back(transVtx.degreesOfFreedom());
               std::vector<float> iTrigMatchdR(3, 999.0);
               vector<TLorentzVector> TrackTriggMatch;

               for (unsigned int i=0; i < iTwoMuTr.size(); i++) {
                  if(i <2 ){	
                     reco::MuonRef tmpRef(muons, iTwoMuTr.at(i));
                     double tmpRef_p = tmpRef->p();
                     TLorentzVector tmpP4(tmpRef->px(), tmpRef->py(), tmpRef->pz(),
                                          sqrt(pow(tmpRef_p,2)+pow(PDGInfo::mu_mass(),2)));
                     TrackTriggMatch.push_back(tmpP4);
                  } else {
                     TrackRef tmpRef(trackCollection, iTwoMuTr.at(i));
                     double tmpRef_p = tmpRef->p();
                     TLorentzVector tmpP4(tmpRef->px(), tmpRef->py(), tmpRef->pz(),
                                          sqrt(pow(tmpRef_p,2)+pow(PDGInfo::mu_mass(),2))); // use muon hypothesis for the track
                     TrackTriggMatch.push_back(tmpP4);
                     dump_track_index_to_fill.push_back(iTwoMuTr.at(i));
                  }
               }
               TriggerMatch(triggerSummary,  TrackTriggMatch , TriggerMuonMatchingdr_, iTrigMatchdR);
               TwoMuonsTrack_TriggerMatch_dR.push_back(iTrigMatchdR);
            }
         }
      }
   }
   return TwoMuonsTrack_idx.size();
}

int T3MNtuple::fillTwoMuonsAndTracks(const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<vector<pat::Muon> >& muons,
      const Handle<TrackCollection>& tracks,
      const Handle<BeamSpot>& beamSpotHandle,
      const Handle<vector<pat::TriggerObjectStandAlone>>& triggerObjects,
      const TriggerNames& triggerNames)
{


   std::vector<std::vector<unsigned int> > PreselectedTwoMuonsTrackCollection  = findTwoMuonsAndTrackCandidates(iEvent, iSetup, beamSpotHandle, muons, trackCollection);
   if(DEBUG)  std::cout<<" PreselectedTwoMuonsTrackCollection  "<< PreselectedTwoMuonsTrackCollection.size()<< std::endl;
   if(PreselectedTwoMuonsTrackCollection.size()==0){
      return 0;          //No two muons + track candidate found!
   }



   for ( auto &iTwoMuTr :  PreselectedTwoMuonsTrackCollection ) {
      vector<TransientTrack> t_trks;
      TransientVertex transVtx;
      ESHandle<TransientTrackBuilder> theB;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
      pat::MuonRef Muon1(muons, iTwoMuTr.at(0));
      pat::MuonRef Muon2(muons, iTwoMuTr.at(1));

      TLorentzVector mv1,mv2;

      mv1.SetPtEtaPhiM(Muon1->pt(), Muon1->eta(), Muon1->phi(), 0.106);
      mv2.SetPtEtaPhiM(Muon2->pt(), Muon2->eta(), Muon2->phi(), 0.106);

      TrackRef track1 = Muon1->innerTrack();
      TrackRef track2 = Muon2->innerTrack();
      TrackRef track3 = TrackRef(trackCollection, iTwoMuTr.at(2));

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



      if(FitOk){
         if(transVtx.totalChiSquared() < 100){
            if( (mv1 + mv2).M() < phimassmin_)// || (mv1 + mv2).M() > phimassmax_)   // renmove that
            {
               TwoMuonsTrack_idx.push_back(iTwoMuTr);
               TwoMuonsTrack_SV_Chi2.push_back(transVtx.totalChiSquared());
               TwoMuonsTrack_SV_NDF.push_back(transVtx.degreesOfFreedom());
               std::vector<float> iTrigMatchdR(3, 999.0);
               vector<TLorentzVector> TrackTriggMatch;

               for (unsigned int i=0; i < iTwoMuTr.size(); i++) {
                  if(i <2 ){	
                     pat::MuonRef tmpRef(muons, iTwoMuTr.at(i));
                     double tmpRef_p = tmpRef->p();
                     TLorentzVector tmpP4(tmpRef->px(), tmpRef->py(), tmpRef->pz(),
                                          sqrt(pow(tmpRef_p,2)+pow(PDGInfo::mu_mass(),2)));
                     TrackTriggMatch.push_back(tmpP4);
                  } else {
                     TrackRef tmpRef(trackCollection, iTwoMuTr.at(i));
                     double tmpRef_p = tmpRef->p();
                     TLorentzVector tmpP4(tmpRef->px(), tmpRef->py(), tmpRef->pz(),
                                          sqrt(pow(tmpRef_p,2)+pow(PDGInfo::mu_mass(),2))); // use muon hypothesis for the track
                     TrackTriggMatch.push_back(tmpP4);
                     dump_track_index_to_fill.push_back(iTwoMuTr.at(i));
                  }
               }
               TriggerMatch(iEvent, triggerObjects,  triggerNames, TrackTriggMatch , TriggerMuonMatchingdr_, iTrigMatchdR);
               TwoMuonsTrack_TriggerMatch_dR.push_back(iTrigMatchdR);
            }
         }
      }
   }
   return TwoMuonsTrack_idx.size();
}
