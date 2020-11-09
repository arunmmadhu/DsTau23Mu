#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"

   std::vector<std::vector<unsigned int> > 
      T3MNtuple::findThreeMuonsCandidates(const edm::Event& iEvent,
                                          const edm::EventSetup& iSetup,
                                          const Handle<BeamSpot>& beamSpotHandle,
                                          const Handle<vector<reco::Muon> >& muonCollection,
                                          const Handle<TrackCollection>& trackCollection)
{

   BeamSpot bs;
   bs = *beamSpotHandle;

   int Muon_index = 0;
   std::vector<unsigned int> preselected_muon_idx;
   std::vector<std::vector<unsigned int> > ThreeMuonsCollection;
   for (vector<reco::Muon>::const_iterator iMuon = muonCollection->begin(); iMuon != muonCollection->end(); ++iMuon, Muon_index++) {
      reco::MuonRef RefMuon(muonCollection, Muon_index);
      if(AcceptedMuon(RefMuon))preselected_muon_idx.push_back(Muon_index);
   }

   if(preselected_muon_idx.size() > 2){
      for(size_t i = 0; i < preselected_muon_idx.size()-1; ++ i){

         std::vector<unsigned int> dump_index;
         reco::MuonRef  Muon1(muonCollection, preselected_muon_idx.at(i));
         for(size_t j = i+1; j < preselected_muon_idx.size(); ++ j){
            MuonRef  Muon2(muonCollection, preselected_muon_idx.at(j));
            double dz_12 = abs(Muon2->vz()-Muon1->vz());  //  INFN
            double dr_12 = deltaR(Muon1->eta(), Muon1->phi(), Muon2->eta(), Muon2->phi());                                                //   not far from each other
            if(dz_12>0.5 ||  dr_12>0.8)continue; // - to be checked  -  this is previsou req.

            if(j<preselected_muon_idx.size()-1){
               for(size_t k = j+1; k < preselected_muon_idx.size(); ++ k){
                  reco::MuonRef  Muon3(muonCollection, preselected_muon_idx.at(k));
                  size_t number_of_muons_pt2p5 = 0;
                  if(Muon1->pt()>2.5)number_of_muons_pt2p5++;
                  if(Muon2->pt()>2.5)number_of_muons_pt2p5++;
                  if(Muon3->pt()>2.5)number_of_muons_pt2p5++;

                  if(Muon1->pt() < 1 or Muon2->pt() <1 or Muon3->pt()<1)std::cout<<"Wrong pt!!!"<< std::endl;

                  //	    if(number_of_muons_pt2p5<2)continue; 
                  //	    double dz_23 = abs(Muon3->innerTrack()->dz(beamSpotHandle->position())-Muon2->innerTrack()->dz(beamSpotHandle->position()));
                  //	    double dz_31 = abs(Muon3->innerTrack()->dz(beamSpotHandle->position())-Muon1->innerTrack()->dz(beamSpotHandle->position()));

                  double dz_23 = abs(Muon3->vz() - Muon2->vz());
                  double dz_31 = abs(Muon3->vz() - Muon1->vz());

                  double dr_23 = deltaR(Muon3->eta(), Muon3->phi(), Muon2->eta(), Muon2->phi());
                  double dr_31 = deltaR(Muon3->eta(), Muon3->phi(), Muon1->eta(), Muon1->phi());

                  if(dr_23>0.8 || dr_31>0.8)continue; 
                  if(dz_23>0.5 || dz_31>0.5)continue; 
                  if(abs(Muon1->charge()+Muon2->charge()+Muon3->charge())>1.1)continue;
                  dump_index.push_back(preselected_muon_idx.at(i));
                  dump_index.push_back(preselected_muon_idx.at(j));
                  dump_index.push_back(preselected_muon_idx.at(k));
                  ThreeMuonsCollection.push_back(dump_index);
                  dump_index.clear();

               }
            }
         }
      }
   }

   return ThreeMuonsCollection;
}

   std::vector<std::vector<unsigned int> > 
      T3MNtuple::findThreeMuonsCandidates(const edm::Event& iEvent,
                                          const edm::EventSetup& iSetup,
                                          const Handle<BeamSpot>& beamSpotHandle,
                                          const Handle<vector<pat::Muon> >& muonCollection,
                                          const Handle<TrackCollection>& trackCollection)
{

   BeamSpot bs;
   bs = *beamSpotHandle;

   int Muon_index = 0;
   std::vector<unsigned int> preselected_muon_idx;
   std::vector<std::vector<unsigned int> > ThreeMuonsCollection;
   for (vector<pat::Muon>::const_iterator iMuon = muonCollection->begin(); iMuon != muonCollection->end(); ++iMuon, Muon_index++) {
      pat::MuonRef RefMuon(muonCollection, Muon_index);
      if(AcceptedMuon(RefMuon))preselected_muon_idx.push_back(Muon_index);
   }

   if(preselected_muon_idx.size() > 2){
      for(size_t i = 0; i < preselected_muon_idx.size()-1; ++ i){

         std::vector<unsigned int> dump_index;
         pat::MuonRef  Muon1(muonCollection, preselected_muon_idx.at(i));
         for(size_t j = i+1; j < preselected_muon_idx.size(); ++ j){
            pat::MuonRef  Muon2(muonCollection, preselected_muon_idx.at(j));
            double dz_12 = abs(Muon2->vz()-Muon1->vz());  //  INFN
            double dr_12 = deltaR(Muon1->eta(), Muon1->phi(), Muon2->eta(), Muon2->phi());                                                //   not far from each other
            if(dz_12>0.5 ||  dr_12>0.8)continue; // - to be checked  -  this is previsou req.

            if(j<preselected_muon_idx.size()-1){
               for(size_t k = j+1; k < preselected_muon_idx.size(); ++ k){
                  pat::MuonRef  Muon3(muonCollection, preselected_muon_idx.at(k));
                  size_t number_of_muons_pt2p5 = 0;
                  if(Muon1->pt()>2.5)number_of_muons_pt2p5++;
                  if(Muon2->pt()>2.5)number_of_muons_pt2p5++;
                  if(Muon3->pt()>2.5)number_of_muons_pt2p5++;

                  if(Muon1->pt() < 1 or Muon2->pt() <1 or Muon3->pt()<1)std::cout<<"Wrong pt!!!"<< std::endl;

                  //	    if(number_of_muons_pt2p5<2)continue; 
                  //	    double dz_23 = abs(Muon3->innerTrack()->dz(beamSpotHandle->position())-Muon2->innerTrack()->dz(beamSpotHandle->position()));
                  //	    double dz_31 = abs(Muon3->innerTrack()->dz(beamSpotHandle->position())-Muon1->innerTrack()->dz(beamSpotHandle->position()));

                  double dz_23 = abs(Muon3->vz() - Muon2->vz());
                  double dz_31 = abs(Muon3->vz() - Muon1->vz());

                  double dr_23 = deltaR(Muon3->eta(), Muon3->phi(), Muon2->eta(), Muon2->phi());
                  double dr_31 = deltaR(Muon3->eta(), Muon3->phi(), Muon1->eta(), Muon1->phi());

                  if(dr_23>0.8 || dr_31>0.8)continue; 
                  if(dz_23>0.5 || dz_31>0.5)continue; 
                  if(abs(Muon1->charge()+Muon2->charge()+Muon3->charge())>1.1)continue;
                  dump_index.push_back(preselected_muon_idx.at(i));
                  dump_index.push_back(preselected_muon_idx.at(j));
                  dump_index.push_back(preselected_muon_idx.at(k));
                  ThreeMuonsCollection.push_back(dump_index);
                  dump_index.clear();

               }
            }
         }
      }
   }

   return ThreeMuonsCollection;
}
