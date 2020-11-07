#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"

void T3MNtuple::fillTracks( const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<TrackCollection>& trackCollection)
{

   unsigned int Track_index(0);
   unsigned int sel_track_index = 0;

   TwoMuonsTrack_Trackindex.resize(TwoMuonsTrack_idx.size());

   std::vector<reco::Track>::const_iterator trIt  = trackCollection->begin();
   std::vector<reco::Track>::const_iterator trEnd = trackCollection->end();

   for (; trIt != trEnd; ++trIt, Track_index++) 
   {
      if(find(dump_track_index_to_fill.begin(), dump_track_index_to_fill.end(), Track_index) !=  dump_track_index_to_fill.end())
      {
         std::vector<double> iTrack_p4;
         std::vector<double> iTrack_poca;
         const reco::Track track = (*trIt);
         if(isGoodTrack(track)){
            for(unsigned int iTwoMuonsTrack=0;  iTwoMuonsTrack < TwoMuonsTrack_idx.size(); iTwoMuonsTrack++){
               if(find(TwoMuonsTrack_idx.at(iTwoMuonsTrack).begin(), TwoMuonsTrack_idx.at(iTwoMuonsTrack).end(), Track_index) !=  TwoMuonsTrack_idx.at(iTwoMuonsTrack).end()){
                  TwoMuonsTrack_Trackindex.at(iTwoMuonsTrack).push_back(sel_track_index);
               }
            }

            iTrack_p4.push_back(sqrt(pow(track.p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
            iTrack_p4.push_back(track.px());
            iTrack_p4.push_back(track.py());
            iTrack_p4.push_back(track.pz());
            Track_p4.push_back(iTrack_p4);
            Track_normalizedChi2.push_back(track.normalizedChi2());
            Track_numberOfValidHits.push_back(track.numberOfValidHits());
            Track_charge.push_back(track.charge());
            Track_dxy.push_back(track.dxy());
            Track_dz.push_back(track.dz());
            iTrack_poca.push_back(track.vx());
            iTrack_poca.push_back(track.vy());
            iTrack_poca.push_back(track.vz());
            Track_poca.push_back(iTrack_poca);

            Track_dxyError.push_back(track.dxyError());
            Track_dzError.push_back(track.dzError());
            sel_track_index++;
         }
      }
   }
}
