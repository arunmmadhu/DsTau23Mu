#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"

   std::vector<std::vector<unsigned int> > 
T3MNtuple::findTwoMuonsAndTrackCandidates( const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<BeamSpot>& beamSpotHandle,
      const Handle<vector<reco::Muon> >& muons,
      const Handle<TrackCollection>& trackCollection)
{

   BeamSpot bs;
   bs = *beamSpotHandle;

   int Muon_index = 0;

   std::vector<unsigned int> preselected_muon_idx;
   std::vector<std::vector<unsigned int> > TwoMuonsPlusTrackCollection; // note that the track index goes last

   for (vector<reco::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon, Muon_index++) {
      reco::MuonRef RefMuon(muons, Muon_index);
      if(AcceptedMuon(RefMuon)) preselected_muon_idx.push_back(Muon_index);
   }
   if(preselected_muon_idx.size() > 1){
      for(size_t i = 0; i < preselected_muon_idx.size()-1; ++ i){
         std::vector<unsigned int> dump_index;
         reco::MuonRef  Muon1(muons, preselected_muon_idx.at(i));
         for(size_t j = i+1; j < preselected_muon_idx.size(); ++ j){
            reco::MuonRef  Muon2(muons, preselected_muon_idx.at(j));
            double dz_12 = abs(Muon2->vz()-Muon1->vz());  // like NFN for sync
            // double dz_12 = abs(Muon2->innerTrack()->dz(beamSpotHandle->position())-Muon1->innerTrack()->dz(beamSpotHandle->position()));  //   Check that two muons are
            double dr_12 = deltaR(Muon1->eta(), Muon1->phi(), Muon2->eta(), Muon2->phi());                                                //   not far from each othe
            // if(abs(Muon1->charge() + Muon2->charge()) !=0 ) continue;

            if(dz_12<0.5 &&  dr_12<0.8)
            { // - to be checked
               unsigned int Track_index = 0;
               for (reco::TrackCollection::const_iterator iTrack = trackCollection->begin(); iTrack != trackCollection->end(); ++iTrack, Track_index++)
               {
                  const reco::Track track = (*iTrack);
                  if(isGoodTrack(track))
                  {
                     double dz23 = fabs(track.vz()  - Muon2->vz());  // like INFN
                     double dz31 = fabs(track.vz()  - Muon1->vz());  // like INFN

                     double dr23 = deltaR(track.eta(), track.phi(), Muon2->eta(), Muon2->phi());
                     double dr31 = deltaR(track.eta(), track.phi(), Muon1->eta(), Muon1->phi());

                     if(dr23 > 0.8  || dr31 > 0.8 )  continue;
                     if(dr23 < 0.01 || dr31 < 0.01)  continue;
                     if(dz23 > 0.5  || dz31 > 0.5 )  continue;

                     if( abs(Muon1->charge() + Muon2->charge() + track.charge())>1.1 ) continue;  // check the charge

                     dump_index.push_back(preselected_muon_idx.at(i));
                     dump_index.push_back(preselected_muon_idx.at(j));
                     dump_index.push_back(Track_index);
                     TwoMuonsPlusTrackCollection.push_back(dump_index);
                     dump_index.clear();
                  }
               }
            }
         }
      }
   }
   return TwoMuonsPlusTrackCollection;
}


   std::vector<std::vector<unsigned int> > 
T3MNtuple::findTwoMuonsAndTrackCandidates( const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<BeamSpot>& beamSpotHandle,
      const Handle<vector<pat::Muon> >& muons,
      const Handle<TrackCollection>& trackCollection)
{

   BeamSpot bs;
   bs = *beamSpotHandle;
   
   int Muon_index = 0;

   std::vector<unsigned int> preselected_muon_idx;
   std::vector<std::vector<unsigned int> > TwoMuonsPlusTrackCollection; // note that the track index goes last

   for (vector<pat::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon, Muon_index++) {
      pat::MuonRef RefMuon(muons, Muon_index);
      if(AcceptedMuon(RefMuon)) preselected_muon_idx.push_back(Muon_index);
   }
   if(preselected_muon_idx.size() > 1){
      for(size_t i = 0; i < preselected_muon_idx.size()-1; ++ i){
         std::vector<unsigned int> dump_index;
         pat::MuonRef  Muon1(muons, preselected_muon_idx.at(i));
         for(size_t j = i+1; j < preselected_muon_idx.size(); ++ j){
            pat::MuonRef  Muon2(muons, preselected_muon_idx.at(j));
            double dz_12 = abs(Muon2->vz()-Muon1->vz());  // like NFN for sync
            // double dz_12 = abs(Muon2->innerTrack()->dz(beamSpotHandle->position())-Muon1->innerTrack()->dz(beamSpotHandle->position()));  //   Check that two muons are
            double dr_12 = deltaR(Muon1->eta(), Muon1->phi(), Muon2->eta(), Muon2->phi());                                                //   not far from each othe
            // if(abs(Muon1->charge() + Muon2->charge()) !=0 ) continue;

            if(dz_12<0.5 &&  dr_12<0.8)
            { // - to be checked
               unsigned int Track_index = 0;
               for (reco::TrackCollection::const_iterator iTrack = trackCollection->begin(); iTrack != trackCollection->end(); ++iTrack, Track_index++)
               {
                  const reco::Track track = (*iTrack);
                  if(isGoodTrack(track))
                  {
                     double dz23 = fabs(track.vz()  - Muon2->vz());  // like INFN
                     double dz31 = fabs(track.vz()  - Muon1->vz());  // like INFN

                     double dr23 = deltaR(track.eta(), track.phi(), Muon2->eta(), Muon2->phi());
                     double dr31 = deltaR(track.eta(), track.phi(), Muon1->eta(), Muon1->phi());

                     if(dr23 > 0.8  || dr31 > 0.8 )  continue;
                     if(dr23 < 0.01 || dr31 < 0.01)  continue;
                     if(dz23 > 0.5  || dz31 > 0.5 )  continue;

                     if( abs(Muon1->charge() + Muon2->charge() + track.charge())>1.1 ) continue;  // check the charge

                     dump_index.push_back(preselected_muon_idx.at(i));
                     dump_index.push_back(preselected_muon_idx.at(j));
                     dump_index.push_back(Track_index);
                     TwoMuonsPlusTrackCollection.push_back(dump_index);
                     dump_index.clear();
                  }
               }
            }
         }
      }
   }
   return TwoMuonsPlusTrackCollection;
}
