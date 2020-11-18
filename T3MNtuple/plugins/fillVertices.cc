#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
#include "DsTau23Mu/T3MNtuple/interface/LorentzVectorParticle.h"

void T3MNtuple::fillVertices(const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<vector<reco::Muon> >& muons,
      const Handle<TrackCollection>& trackCollection, 
      const Handle<VertexCollection >& pvs,
      const Handle<std::vector<reco::Vertex> >& svertices,
      const Handle<reco::BeamSpot>& beamSpotHandle)
{

   const BeamSpot& bs = *beamSpotHandle;

   std::vector<std::vector<std::vector<double> > > particles_p4;
   std::vector<std::vector<TransientTrack> > signalTracksCollection;
   ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
   if(ThreeMuons_idx.size()!=0){
      for ( auto &iThreeMuon :  ThreeMuons_idx ) {
         particles_p4.push_back(std::vector<std::vector<double> > ());
         vector<TransientTrack> isignalTracksCollection;
         ESHandle<TransientTrackBuilder> theB;
         iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
         for ( auto &iMuon :  iThreeMuon ) {
            reco::MuonRef Muon(muons, iMuon);

            TrackRef MuonTrack = Muon->innerTrack();
            isignalTracksCollection.push_back(theB->build(MuonTrack));
            std::vector<double> iiparticles_p4;
            iiparticles_p4.push_back(Muon->p4().E());
            iiparticles_p4.push_back(Muon->p4().Px());
            iiparticles_p4.push_back(Muon->p4().Py());
            iiparticles_p4.push_back(Muon->p4().Pz());
            particles_p4.at(particles_p4.size() -1).push_back(iiparticles_p4);
         }
         signalTracksCollection.push_back(isignalTracksCollection);
      }
   }

   if(TwoMuonsTrack_idx.size()!=0){
      std::vector<std::vector<double> > iparticle_p4;
      for ( auto &iTwoMuonsTracks :  TwoMuonsTrack_idx ) {

         vector<TransientTrack> isignalTracksCollection;

         reco::MuonRef Muon1(muons, iTwoMuonsTracks.at(0));
         reco::MuonRef Muon2(muons, iTwoMuonsTracks.at(1));

         TrackRef track1 = Muon1->innerTrack();
         TrackRef track2 = Muon2->innerTrack();
         TrackRef track3 = TrackRef(trackCollection, iTwoMuonsTracks.at(2));

         isignalTracksCollection.push_back(theB->build(track1));
         isignalTracksCollection.push_back(theB->build(track2));
         isignalTracksCollection.push_back(theB->build(track3));
         signalTracksCollection.push_back(isignalTracksCollection);

         std::vector<double> particle1_p4;
         std::vector<double> particle2_p4;
         std::vector<double> particle3_p4;

         particle1_p4.push_back(Muon1->p4().E());       particle2_p4.push_back(Muon2->p4().E());       particle3_p4.push_back(sqrt(pow(track3->p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
         particle1_p4.push_back(Muon1->p4().Px());      particle2_p4.push_back(Muon2->p4().Px());      particle3_p4.push_back(track3->px());
         particle1_p4.push_back(Muon1->p4().Py());      particle2_p4.push_back(Muon2->p4().Py());      particle3_p4.push_back(track3->py());
         particle1_p4.push_back(Muon1->p4().Pz());      particle2_p4.push_back(Muon2->p4().Pz());      particle3_p4.push_back(track3->pz());
         iparticle_p4.push_back(particle1_p4);          iparticle_p4.push_back(particle2_p4);          iparticle_p4.push_back(particle3_p4);
         particles_p4.push_back(iparticle_p4);
      }
   }

   if (DEBUG)   std::cout<<" SignalTracksCollection size:  "<< signalTracksCollection.size() << std::endl;
   unsigned int index(0);
   for ( auto &iTransientTracks :  signalTracksCollection ){
      Vertex_signal_KF_pos.push_back(std::vector<double> ());
      Vertex_signal_KF_cov.push_back(std::vector<double> ());
      Vertex_signal_KF_refittedTracksP4.push_back(std::vector<std::vector<double> >());

      Vertex_signal_AF_pos.push_back(std::vector<double> ());
      ClosestApproachInRPhi cApp12, cApp23, cApp31;
      cApp12.calculate(iTransientTracks[0].initialFreeState(), iTransientTracks[1].initialFreeState());
      cApp23.calculate(iTransientTracks[1].initialFreeState(), iTransientTracks[2].initialFreeState());
      cApp31.calculate(iTransientTracks[2].initialFreeState(), iTransientTracks[0].initialFreeState());
      std::vector<double> iVertex_signal_dca_reco;
      if((cApp12.status()&&cApp23.status()&&cApp31.status())) { 
         iVertex_signal_dca_reco.push_back(cApp12.distance());   // order 12,23,31, max
         iVertex_signal_dca_reco.push_back(cApp23.distance());
         iVertex_signal_dca_reco.push_back(cApp31.distance());
         iVertex_signal_dca_reco.push_back(TMath::Max(dca12_reco, TMath::Max(dca23_reco, dca31_reco)));
      } else {
         iVertex_signal_dca_reco.push_back(-1.);
         iVertex_signal_dca_reco.push_back(-1.);
         iVertex_signal_dca_reco.push_back(-1.);
         iVertex_signal_dca_reco.push_back(-1.);
      }
      Vertex_signal_dca_reco.push_back(iVertex_signal_dca_reco);

      TransientVertex transVtx;
      KalmanVertexFitter kvf(true);
      bool FitOk(true);
      try {
         transVtx = kvf.vertex(iTransientTracks); 
      } catch (...) {
         FitOk = false;
      }
      if (!transVtx.hasRefittedTracks())
         FitOk = false;
      if (transVtx.refittedTracks().size() != iTransientTracks.size())
         FitOk = false;
      TLorentzVector ThreeCandidate(0,0,0,0);

      math::XYZPoint TheSecondaryVertexPoint;
      if(FitOk){
         Vertex_signal_KF_Chi2.push_back(transVtx.totalChiSquared());
         Vertex_signal_KF_pos.at(index).push_back(transVtx.position().x());
         Vertex_signal_KF_pos.at(index).push_back(transVtx.position().y());
         Vertex_signal_KF_pos.at(index).push_back(transVtx.position().z());

         reco::Vertex secondaryVertex = transVtx;
         TMatrixTSym<double> svcov(LorentzVectorParticle::NVertex);
         math::Error<3>::type svCov;
         secondaryVertex.fill(svCov);

         for (int i = 0; i <LorentzVectorParticle::NVertex; i++){
            for (int j = 0; j < LorentzVectorParticle::NVertex; j++) {
               svcov(i, j) = svCov(i, j);
               svcov(j, i) = svCov(i, j);
            }
         }
         for (int i = 0; i < LorentzVectorParticle::NVertex; i++) {
            for (int j = i; j < LorentzVectorParticle::NVertex; j++) {
               Vertex_signal_KF_cov.at(index).push_back(svcov(i, j));
            }
         }

         TheSecondaryVertexPoint = math::XYZPoint(transVtx.position().x(), transVtx.position().y(), transVtx.position().z());
         vector<TransientTrack>::const_iterator trkIt = transVtx.refittedTracks().begin();

         reco::Vertex TripletVtx = reco::Vertex(TheSecondaryVertexPoint, svCov, transVtx.totalChiSquared(), 3, 3);

         VertexState BSstate(bs);
         VertexDistanceXY vertTool2D;
         double BSdistance2D = vertTool2D.distance(BSstate, TripletVtx).value();
         double BSdist_err2D = vertTool2D.distance(BSstate, TripletVtx).error();
         double BSdist_sign2D =vertTool2D.distance(BSstate, TripletVtx).significance();

         Vertex_signal_KF_BS_2Ddistance.push_back(BSdistance2D);
         Vertex_signal_KF_BS_error.push_back(BSdist_err2D);
         Vertex_signal_KF_BS_significance.push_back(BSdist_sign2D);

         for(; trkIt != transVtx.refittedTracks().end(); ++ trkIt) {
            std::vector<double> irefitted_tracks_p4;
            const Track & trkrefit = trkIt->track();
            irefitted_tracks_p4.push_back(sqrt(pow(trkrefit.p(),2.0) + pow(PDGInfo::mu_mass(),2.0)));
            irefitted_tracks_p4.push_back(trkrefit.px());
            irefitted_tracks_p4.push_back(trkrefit.py());
            irefitted_tracks_p4.push_back(trkrefit.pz());
            ThreeCandidate+=TLorentzVector(trkrefit.px(),trkrefit.py(),trkrefit.pz(),sqrt(pow(trkrefit.p(),2.0) + pow(PDGInfo::mu_mass(),2.0)));
            Vertex_signal_KF_refittedTracksP4.at(Vertex_signal_KF_refittedTracksP4.size() -1).push_back(irefitted_tracks_p4);
         }
      }

      bool AFitOk(true);
      AdaptiveVertexFitter avf;
      TransientVertex AdaptivetransVtx;
      avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
      try {
         AdaptivetransVtx = avf.vertex(iTransientTracks);
      } catch (...) {
         AFitOk = false;
      }
      if(!AdaptivetransVtx.isValid()) 
         AFitOk = false;
      if (AFitOk){ 
         Vertex_signal_AF_Chi2.push_back(AdaptivetransVtx.totalChiSquared());
         Vertex_signal_AF_Ndf.push_back(AdaptivetransVtx.degreesOfFreedom());
         Vertex_signal_AF_pos.at(index).push_back(AdaptivetransVtx.position().x());
         Vertex_signal_AF_pos.at(index).push_back(AdaptivetransVtx.position().y());
         Vertex_signal_AF_pos.at(index).push_back(AdaptivetransVtx.position().z());

      } else {
         Vertex_signal_AF_Chi2.push_back(-1);
      }

      //--------------------  Fit each track pair 
      vector<TransientTrack> trackpair12, trackpair23, trackpair31;
      trackpair12.push_back(iTransientTracks.at(0)); trackpair12.push_back(iTransientTracks.at(1));
      trackpair23.push_back(iTransientTracks.at(1)); trackpair23.push_back(iTransientTracks.at(2));
      trackpair31.push_back(iTransientTracks.at(2)); trackpair31.push_back(iTransientTracks.at(0));
      KalmanVertexFitter kvf_trks12(true), kvf_trks23(true), kvf_trks31(true);
      TransientVertex fv_trks12;// = kvf_trks12.vertex(trackpair12);
      TransientVertex fv_trks23;// = kvf_trks23.vertex(trackpair23);
      TransientVertex fv_trks31;// = kvf_trks31.vertex(trackpair31);


      bool Fit1Ok(true);
      try {
         fv_trks12 = kvf_trks12.vertex(trackpair12); 
      } catch (...) {
         Fit1Ok = false;
      }

      bool Fit2Ok(true);
      try {
         fv_trks23 = kvf_trks23.vertex(trackpair23); 
      } catch (...) {
         Fit2Ok = false;
      }

      bool Fit3Ok(true);
      try {
         fv_trks31 = kvf_trks31.vertex(trackpair31); 
      } catch (...) {
         Fit3Ok = false;
      }

      std::vector<double> iVertex_pair_quality;
      std::vector<double> iVertex_pairfit_status;
      iVertex_pairfit_status.push_back(Fit1Ok);
      iVertex_pairfit_status.push_back(Fit2Ok);
      iVertex_pairfit_status.push_back(Fit3Ok);

      std::vector<float> iVertex_Pair12_Pos;
      std::vector<float> iVertex_Pair23_Pos;
      std::vector<float> iVertex_Pair31_Pos;
      if(fv_trks12.isValid())
      {
         iVertex_pair_quality.push_back(fv_trks12.totalChiSquared());
         iVertex_Pair12_Pos.push_back( fv_trks12.position().x());
         iVertex_Pair12_Pos.push_back( fv_trks12.position().y());
         iVertex_Pair12_Pos.push_back( fv_trks12.position().z());
      }

      else
      {
         iVertex_pair_quality.push_back(-1);
         iVertex_Pair12_Pos.push_back( 99.);
         iVertex_Pair12_Pos.push_back( 99.);
         iVertex_Pair12_Pos.push_back( 99.);
      }
      if(fv_trks23.isValid())
      {
         iVertex_pair_quality.push_back(fv_trks23.totalChiSquared());
         iVertex_Pair23_Pos.push_back( fv_trks23.position().x());
         iVertex_Pair23_Pos.push_back( fv_trks23.position().y());
         iVertex_Pair23_Pos.push_back( fv_trks23.position().z());
      }
      else
      {
         iVertex_pair_quality.push_back(-1);
         iVertex_Pair23_Pos.push_back( 99.);
         iVertex_Pair23_Pos.push_back( 99.);
         iVertex_Pair23_Pos.push_back( 99.);
      }
      if(fv_trks31.isValid())
      {
         iVertex_pair_quality.push_back(fv_trks31.totalChiSquared());
         iVertex_Pair31_Pos.push_back( fv_trks31.position().x());
         iVertex_Pair31_Pos.push_back( fv_trks31.position().y());
         iVertex_Pair31_Pos.push_back( fv_trks31.position().z());
      }
      else
      {
         iVertex_pair_quality.push_back(-1);
         iVertex_Pair31_Pos.push_back( 99.);
         iVertex_Pair31_Pos.push_back( 99.);
         iVertex_Pair31_Pos.push_back( 99.);
      }

      Vertex_Pair12_Pos.push_back( iVertex_Pair12_Pos);
      Vertex_Pair23_Pos.push_back( iVertex_Pair23_Pos);
      Vertex_Pair31_Pos.push_back( iVertex_Pair31_Pos);


      Vertex_pair_quality.push_back(iVertex_pair_quality);
      Vertex_pairfit_status.push_back(iVertex_pairfit_status);

      ///////////////////////////////////////////
      //  find here the primary vertex with the best
      //  alignement to the tri-muon 
      double dphi_pv = -1.0;
      unsigned int primaryvertex_index;
      for(unsigned int vertex_index = 0; vertex_index  < pvs->size(); vertex_index++) {
         const Vertex & pvertex = (*pvs)[vertex_index];
         TVector3 Dv3D_reco(transVtx.position().x() - pvertex.x(), transVtx.position().y() - pvertex.y(), transVtx.position().z() - pvertex.z());
         double Cosdphi_3D = Dv3D_reco.Dot(ThreeCandidate.Vect())/(Dv3D_reco.Mag()*ThreeCandidate.Vect().Mag());
         if(Cosdphi_3D>dphi_pv){
            dphi_pv = Cosdphi_3D;
            primaryvertex_index=vertex_index;
         }
      }

      const Vertex & MatchedPrimaryVertex = (*pvs)[primaryvertex_index];
      dump_pv_index_to_fill.push_back(primaryvertex_index);


      std::vector<double>  iprimaryVertex_Pos;
      iprimaryVertex_Pos.push_back(MatchedPrimaryVertex.x());
      iprimaryVertex_Pos.push_back(MatchedPrimaryVertex.y());
      iprimaryVertex_Pos.push_back(MatchedPrimaryVertex.z());
      Vertex_MatchedPrimaryVertex.push_back(iprimaryVertex_Pos);

      double tempdz(99.);
      unsigned int secondbest_primaryvertex_index(0);
      for(unsigned int vertex_index = 0; vertex_index  < pvs->size(); vertex_index++) {
         if(vertex_index == primaryvertex_index) continue;
         const Vertex & temp_pv = (*pvs)[vertex_index];
         if(fabs(temp_pv.z() -  MatchedPrimaryVertex.z()) < tempdz ){
            tempdz = fabs(temp_pv.z() -  MatchedPrimaryVertex.z());
            secondbest_primaryvertex_index = vertex_index;
         }
      }


      const Vertex & SecondBestPrimaryVertex = (*pvs)[secondbest_primaryvertex_index];
      std::vector<double> iSecondBestprimaryVertex_Pos;
      if( pvs->size()>1){
         iSecondBestprimaryVertex_Pos.push_back(SecondBestPrimaryVertex.x());
         iSecondBestprimaryVertex_Pos.push_back(SecondBestPrimaryVertex.y());
         iSecondBestprimaryVertex_Pos.push_back(SecondBestPrimaryVertex.z());
      }else{
         iSecondBestprimaryVertex_Pos.push_back(-99.);
         iSecondBestprimaryVertex_Pos.push_back(-99.);
         iSecondBestprimaryVertex_Pos.push_back(-99.);
      }
      Vertex_SecondBestPrimaryVertex.push_back(iSecondBestprimaryVertex_Pos);


      int NParticlesComingFromPV(0);
      vector<TransientTrack> primaryvertexTransientTracks;// remove muon candidates from the PV to perform refit
      for(Vertex::trackRef_iterator itk = MatchedPrimaryVertex.tracks_begin(); itk != MatchedPrimaryVertex.tracks_end(); itk++) {

         if((**itk).pt()>1) {
            if(deltaR(iTransientTracks.at(0).track().eta(), iTransientTracks.at(0).track().phi(), (**itk).eta(), (**itk).phi())<0.01){NParticlesComingFromPV++;continue;}
            if(deltaR(iTransientTracks.at(1).track().eta(), iTransientTracks.at(1).track().phi(), (**itk).eta(), (**itk).phi())<0.01){NParticlesComingFromPV++;continue;}
            if(deltaR(iTransientTracks.at(2).track().eta(), iTransientTracks.at(2).track().phi(), (**itk).eta(), (**itk).phi())<0.01){NParticlesComingFromPV++;continue;}
         }

         primaryvertexTransientTracks.push_back(theB->build(**itk));

         std::vector<float> iIsolationBranch_Track_p4;

         iIsolationBranch_Track_p4.push_back(sqrt(pow((**itk).p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
         iIsolationBranch_Track_p4.push_back((**itk).px());
         iIsolationBranch_Track_p4.push_back((**itk).py());
         iIsolationBranch_Track_p4.push_back((**itk).pz());

         IsolationBranch_Trackp4.at(IsolationBranch_Trackp4.size() - 1).push_back(iIsolationBranch_Track_p4);
      }

      Vertex_NMuonsAssocWithPV.push_back(NParticlesComingFromPV);
      KalmanVertexFitter pv_fit(true);
      bool FitPVOk(true);
      TransientVertex pvvertex;

      if(primaryvertexTransientTracks.size() >1){
         try {
            pvvertex =  pv_fit.vertex(primaryvertexTransientTracks);
         } catch (...) {
            FitPVOk = false;
         }
      }

      Vertex_RefitPVisValid.push_back(pvvertex.isValid());
      std::vector<double> iRefitprimaryVertex_Pos;
      if(FitPVOk && pvvertex.isValid()){
         iRefitprimaryVertex_Pos.push_back(pvvertex.position().x());
         iRefitprimaryVertex_Pos.push_back(pvvertex.position().y());
         iRefitprimaryVertex_Pos.push_back(pvvertex.position().z());
      }
      Vertex_MatchedRefitPrimaryVertex.push_back(iRefitprimaryVertex_Pos);


      TMatrixTSym<double> pvcov(LorentzVectorParticle::NVertex);

      math::Error<3>::type pvCov;
      MatchedPrimaryVertex.fill(pvCov);

      for (int i = 0; i <LorentzVectorParticle::NVertex; i++){
         for (int j = 0; j < LorentzVectorParticle::NVertex; j++) {
            pvcov(i, j) = pvCov(i, j);
            pvcov(j, i) = pvCov(i, j);
         }
      }

      std::vector<double>  pv_cov;     
      for (int i = 0; i < LorentzVectorParticle::NVertex; i++) {
         for (int j = i; j < LorentzVectorParticle::NVertex; j++) {
            pv_cov.push_back(pvcov(i, j));
         }
      }

      Vertex_MatchedRefitPrimaryVertex_covariance.push_back(pv_cov);

      Vertex final_pv = MatchedPrimaryVertex;  
      if(pvvertex.isValid()) final_pv = Vertex(pvvertex);

      math::XYZPoint pvPoint = math::XYZPoint(final_pv.x(), final_pv.y(), final_pv.z());
      math::XYZPoint bsPoint = math::XYZPoint(beamSpotHandle->position().x(), beamSpotHandle->position().y(), beamSpotHandle->position().z());



      std::vector<double> iVertex_d0BeamSpot_reco;
      iVertex_d0BeamSpot_reco.push_back(abs(iTransientTracks.at(0).track().dxy(bsPoint)));
      iVertex_d0BeamSpot_reco.push_back(abs(iTransientTracks.at(1).track().dxy(bsPoint)));
      iVertex_d0BeamSpot_reco.push_back(abs(iTransientTracks.at(2).track().dxy(bsPoint)));
      Vertex_d0BeamSpot_reco.push_back(iVertex_d0BeamSpot_reco);


      std::vector<double> iVertex_d0BeamSpot_reco_sig;
      double d0ErrorToBs_1  = sqrt( iTransientTracks.at(0).track().d0Error() * iTransientTracks.at(0).track().d0Error() +
            0.5*  beamSpotHandle->BeamWidthX()* beamSpotHandle->BeamWidthX()+
            0.5*  beamSpotHandle->BeamWidthY()* beamSpotHandle->BeamWidthY() );

      double d0ErrorToBs_2  = sqrt( iTransientTracks.at(1).track().d0Error() * iTransientTracks.at(1).track().d0Error() +
            0.5*  beamSpotHandle->BeamWidthX()* beamSpotHandle->BeamWidthX()+
            0.5*  beamSpotHandle->BeamWidthY()* beamSpotHandle->BeamWidthY() );

      double d0ErrorToBs_3  = sqrt( iTransientTracks.at(2).track().d0Error() * iTransientTracks.at(2).track().d0Error() +
            0.5*  beamSpotHandle->BeamWidthX()* beamSpotHandle->BeamWidthX()+
            0.5*  beamSpotHandle->BeamWidthY()* beamSpotHandle->BeamWidthY() );



      if(d0ErrorToBs_1!=0){  iVertex_d0BeamSpot_reco_sig.push_back( abs(iTransientTracks.at(0).track().dxy(bsPoint)) / d0ErrorToBs_1);} else {iVertex_d0BeamSpot_reco_sig.push_back(-1);}
      if(d0ErrorToBs_2!=0){  iVertex_d0BeamSpot_reco_sig.push_back( abs(iTransientTracks.at(1).track().dxy(bsPoint)) / d0ErrorToBs_2);} else {iVertex_d0BeamSpot_reco_sig.push_back(-1);}
      if(d0ErrorToBs_3!=0){  iVertex_d0BeamSpot_reco_sig.push_back( abs(iTransientTracks.at(1).track().dxy(bsPoint)) / d0ErrorToBs_3);} else {iVertex_d0BeamSpot_reco_sig.push_back(-1);}

      Vertex_d0BeamSpot_reco_sig.push_back(iVertex_d0BeamSpot_reco_sig);


      std::vector<double> iVertex_d0SV_reco;
      iVertex_d0SV_reco.push_back(abs(iTransientTracks.at(0).track().dxy(TheSecondaryVertexPoint)));
      iVertex_d0SV_reco.push_back(abs(iTransientTracks.at(1).track().dxy(TheSecondaryVertexPoint)));
      iVertex_d0SV_reco.push_back(abs(iTransientTracks.at(2).track().dxy(TheSecondaryVertexPoint)));
      Vertex_d0SV_reco.push_back(iVertex_d0SV_reco);


      std::vector<double> iVertex_dzSV_reco;
      iVertex_dzSV_reco.push_back(abs(iTransientTracks.at(0).track().dz(TheSecondaryVertexPoint)));
      iVertex_dzSV_reco.push_back(abs(iTransientTracks.at(1).track().dz(TheSecondaryVertexPoint)));
      iVertex_dzSV_reco.push_back(abs(iTransientTracks.at(2).track().dz(TheSecondaryVertexPoint)));
      Vertex_dzSV_reco.push_back(iVertex_dzSV_reco);



      std::vector<double> iVertex_d0_reco;
      iVertex_d0_reco.push_back(abs(iTransientTracks.at(0).track().dxy(pvPoint)));
      iVertex_d0_reco.push_back(abs(iTransientTracks.at(1).track().dxy(pvPoint)));
      iVertex_d0_reco.push_back(abs(iTransientTracks.at(2).track().dxy(pvPoint)));
      Vertex_d0_reco.push_back(iVertex_d0_reco);


      std::vector<double> iVertex_dz_reco;
      iVertex_dz_reco.push_back(abs(iTransientTracks.at(0).track().dz(pvPoint)));
      iVertex_dz_reco.push_back(abs(iTransientTracks.at(1).track().dz(pvPoint)));
      iVertex_dz_reco.push_back(abs(iTransientTracks.at(2).track().dz(pvPoint)));
      Vertex_dz_reco.push_back(iVertex_dz_reco);



      TLorentzVector LV1=TLorentzVector(particles_p4.at(index).at(0).at(1), particles_p4.at(index).at(0).at(2), particles_p4.at(index).at(0).at(3),particles_p4.at(index).at(0).at(0));
      TLorentzVector LV2=TLorentzVector(particles_p4.at(index).at(1).at(1), particles_p4.at(index).at(1).at(2), particles_p4.at(index).at(1).at(3),particles_p4.at(index).at(1).at(0));
      TLorentzVector LV3=TLorentzVector(particles_p4.at(index).at(2).at(1), particles_p4.at(index).at(2).at(2), particles_p4.at(index).at(2).at(3),particles_p4.at(index).at(2).at(0));
      TLorentzVector LVTau = LV1 + LV2 + LV3;
      m3mu_reco = LVTau.M();

      GlobalVector dir1(particles_p4.at(index).at(0).at(1), particles_p4.at(index).at(0).at(2), particles_p4.at(index).at(0).at(3));
      GlobalVector dir2(particles_p4.at(index).at(1).at(1), particles_p4.at(index).at(1).at(2), particles_p4.at(index).at(1).at(3));
      GlobalVector dir3(particles_p4.at(index).at(2).at(1), particles_p4.at(index).at(2).at(2), particles_p4.at(index).at(2).at(3));

      std::pair<bool, Measurement1D> ip2d_1 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(0), dir1, final_pv);
      std::pair<bool, Measurement1D> ip2d_2 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(1), dir2, final_pv);
      std::pair<bool, Measurement1D> ip2d_3 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(2), dir3, final_pv);


      std::pair<bool, Measurement1D> ip2dSV_1 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(0), dir1, transVtx);
      std::pair<bool, Measurement1D> ip2dSV_2 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(1), dir2, transVtx);
      std::pair<bool, Measurement1D> ip2dSV_3 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(2), dir3, transVtx);


      //  std::pair<bool, Measurement1D> ipBS2d_1 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(0), dir1, );
      //  std::pair<bool, Measurement1D> ipBS2d_2 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(1), dir2, final_pv);
      //  std::pair<bool, Measurement1D> ipBS2d_3 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(2), dir3, final_pv);

      std::vector<double> iVertex_d0sig_reco;
      if(ip2d_1.first){ iVertex_d0sig_reco.push_back(abs(ip2d_1.second.value()/ip2d_1.second.error()));} else{iVertex_d0sig_reco.push_back(-1);}
      if(ip2d_2.first){ iVertex_d0sig_reco.push_back(abs(ip2d_2.second.value()/ip2d_2.second.error()));} else{iVertex_d0sig_reco.push_back(-1);}
      if(ip2d_3.first){ iVertex_d0sig_reco.push_back(abs(ip2d_3.second.value()/ip2d_3.second.error()));} else{iVertex_d0sig_reco.push_back(-1);}
      Vertex_d0sig_reco.push_back(iVertex_d0sig_reco);


      std::vector<double> iVertex_d0sigSV_reco;
      if(ip2dSV_1.first){ iVertex_d0sigSV_reco.push_back(abs(ip2dSV_1.second.value()/ip2dSV_1.second.error()));} else{iVertex_d0sigSV_reco.push_back(-1);}
      if(ip2dSV_2.first){ iVertex_d0sigSV_reco.push_back(abs(ip2dSV_2.second.value()/ip2dSV_2.second.error()));} else{iVertex_d0sigSV_reco.push_back(-1);}
      if(ip2dSV_3.first){ iVertex_d0sigSV_reco.push_back(abs(ip2dSV_3.second.value()/ip2dSV_3.second.error()));} else{iVertex_d0sigSV_reco.push_back(-1);}
      Vertex_d0sigSV_reco.push_back(iVertex_d0sigSV_reco);


      TVector3 dv2D_reco(-final_pv.position().x() + TheSecondaryVertexPoint.x(), -final_pv.position().y() + TheSecondaryVertexPoint.y(), 0);
      TVector3 vtauxy(ThreeCandidate.Px(), ThreeCandidate.Py(), 0);
      fv_cosdphi = dv2D_reco.Dot(vtauxy)/(dv2D_reco.Perp()*vtauxy.Perp());
      VertexDistanceXY vdistXY;
      Measurement1D distXY = vdistXY.distance(Vertex(transVtx), final_pv);
      std::vector<double> iVertex_2Ddisplacement;
      iVertex_2Ddisplacement.push_back(distXY.value());
      iVertex_2Ddisplacement.push_back(distXY.significance());
      iVertex_2Ddisplacement.push_back(distXY.value()*fv_cosdphi * m3mu_reco/vtauxy.Perp());
      Vertex_2Ddisplacement.push_back(iVertex_2Ddisplacement);


      TVector3 vtauxyz(ThreeCandidate.Px(), ThreeCandidate.Py(), ThreeCandidate.Pz());
      TVector3 dv3D_reco(-final_pv.position().x() + TheSecondaryVertexPoint.x(), -final_pv.position().y() + TheSecondaryVertexPoint.y(), -final_pv.position().z() + TheSecondaryVertexPoint.z());
      fv_cosdphi3D = dv3D_reco.Dot(vtauxyz)/(dv3D_reco.Mag()*vtauxyz.Mag());
      VertexDistance3D dist;

      std::vector<double> iVertex_3Ddisplacement;
      iVertex_3Ddisplacement.push_back(dist.distance(Vertex(transVtx), final_pv).value());
      iVertex_3Ddisplacement.push_back(dist.distance(Vertex(transVtx), final_pv).significance());
      iVertex_3Ddisplacement.push_back(fv_d3D*fv_cosdphi3D*m3mu_reco/ThreeCandidate.P());
      Vertex_3Ddisplacement.push_back(iVertex_3Ddisplacement);


      ///////////////////////////////////////////////////////////////
      //    Here fill the isolation

      //  former Isolation
      float minmuon_pt(999.), maxmuon_dr(0.);

      if(LV1.Pt()<minmuon_pt)minmuon_pt=LV1.Pt();
      if(LV2.Pt()<minmuon_pt)minmuon_pt=LV2.Pt();
      if(LV3.Pt()<minmuon_pt)minmuon_pt=LV3.Pt();

      float  drLV1Tau(deltaR(LV1.Eta(), LV1.Phi(), LVTau.Eta(), LVTau.Phi()));
      float  drLV2Tau(deltaR(LV2.Eta(), LV2.Phi(), LVTau.Eta(), LVTau.Phi()));
      float  drLV3Tau(deltaR(LV3.Eta(), LV3.Phi(), LVTau.Eta(), LVTau.Phi()));

      maxmuon_dr = TMath::Max(drLV1Tau, TMath::Max(drLV2Tau,drLV3Tau));

      float sumptalltracks(0.), sumalltracks(0.), mindist(99.);
      float sumptalltracks05(0.), sumalltracks05(0.), mindist05(99.), sumalltracks_b(0.);
      float pt_trk_1(0.), pt_trk_2(0.), pt_trk_3(0.);
      float N_trk_1(0.),N_trk_2(0.),N_trk_3(0.), N_trk_total(0.);
      float N_trk0p1(0.), N_trk0p2(0.), N_trk0p5(0.), maxdxy(0.);
      float relative_iso(0.), relative_iso05(0.), relative_mu1_iso(0.),relative_mu2_iso(0.),relative_mu3_iso(0.);
      float relative_maxiso(0.);


      //--------------------------- Isolation Branch

      std::vector<TLorentzVector>  MuLVs;
      MuLVs.push_back(LV1);
      MuLVs.push_back(LV2);
      MuLVs.push_back(LV3);


      std::vector<int> sortedindices = SortByPt(MuLVs);


      IsolationTrack_p4.push_back(std::vector<std::vector<float> >());
      IsolationTrack_VertexWithSignalMuonIsValid.push_back(std::vector<std::vector<int> >());
      IsolationTrack_VertexWithSignalMuonChi2.push_back(std::vector<std::vector<float> >());
      IsolationTrack_VertexWithSignalMuonPosition.push_back(std::vector<std::vector<float> >());

      std::vector<float> iIsolationTrack_dxySV;
      std::vector<float> iIsolationTrack_dzSV;
      std::vector<float> iIsolationTrack_dxyPV;
      std::vector<float> iIsolationTrack_dzPV;
      std::vector<float> iIsolationTrack_DocaMu1;
      std::vector<float> iIsolationTrack_DocaMu2;
      std::vector<float> iIsolationTrack_DocaMu3;
      std::vector<int>   iIsolationTrack_charge;
      std::vector<int>   iIsolationTrack_isHighPurity;
      std::vector<int>   IsolationTracksIndices;

      vector<TransientTrack> iIsolationTracksCollection;
      ESHandle<TransientTrackBuilder> TheB;
      std::vector<std::vector<int> >  VertexWithSignalMuonIsValid;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TheB);

      for(size_t i = 0; i < trackCollection->size(); i++) {
         const Track & t = (*trackCollection)[i];


         if(!(t.quality(TrackBase::tight)))continue; //-- this might be weaker
         if(deltaR(LV1.Eta(), LV1.Phi(), t.eta(), t.phi())<0.01)continue;
         if(deltaR(LV2.Eta(), LV2.Phi(), t.eta(), t.phi())<0.01)continue;
         if(deltaR(LV3.Eta(), LV3.Phi(), t.eta(), t.phi())<0.01)continue;
         if(index>=ThreeMuons_idx.size()) break;  // <----- fillinf isolation vertexing only for a signal candidate
         //if(abs(t.dz(pvPoint))< 0.5 && t.quality(TrackBase::tight) && sqrt(t.px()*t.px() + t.py()*t.py() ) > 0.5){//  && deltaR(t.eta(), t.phi(), LVTau.Eta(), LVTau.Phi()) < 1.){}}
         if(abs(t.dz(pvPoint))< 0.5 && sqrt(t.px()*t.px() + t.py()*t.py() ) > 0.4  && deltaR(t.eta(), t.phi(), LVTau.Eta(), LVTau.Phi()) < 1.){

            std::vector<int>   iIsolationTrack_VertexWithSignalMuonIsValid;
            std::vector<float> iIsolationTrack_VertexWithSignalMuonChi2;
            std::vector<float> iIsolationTrack_VertexWithSignalMuonPosition;
            IsolationTracksIndices.push_back(i);
            std::vector<float> iIsolation_Track_p4;
            iIsolation_Track_p4.push_back(sqrt(pow(t.p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
            iIsolation_Track_p4.push_back(t.px());
            iIsolation_Track_p4.push_back(t.py());
            iIsolation_Track_p4.push_back(t.pz());

            iIsolationTrack_charge.push_back(t.charge());
            iIsolationTrack_isHighPurity.push_back(t.quality(TrackBase::highPurity));

            IsolationTrack_p4.at(IsolationTrack_p4.size() - 1).push_back(iIsolation_Track_p4);

            iIsolationTrack_dxySV.push_back(t.dxy(TheSecondaryVertexPoint));
            iIsolationTrack_dzSV.push_back(t.dz(TheSecondaryVertexPoint));

            iIsolationTrack_dxyPV.push_back(t.dxy(pvPoint));
            iIsolationTrack_dzPV.push_back(t.dz(pvPoint));


            ClosestApproachInRPhi DocaMuon1, DocaMuon2, DocaMuon3;

            DocaMuon1.calculate(theB->build(t).initialFreeState(), iTransientTracks.at(sortedindices.at(0)).initialFreeState());
            DocaMuon2.calculate(theB->build(t).initialFreeState(), iTransientTracks.at(sortedindices.at(1)).initialFreeState());
            DocaMuon3.calculate(theB->build(t).initialFreeState(), iTransientTracks.at(sortedindices.at(2)).initialFreeState());


            if(DocaMuon1.status()){ iIsolationTrack_DocaMu1.push_back(DocaMuon1.distance());}
            else iIsolationTrack_DocaMu1.push_back(-1);

            if(DocaMuon2.status()){ iIsolationTrack_DocaMu2.push_back(DocaMuon2.distance());}
            else iIsolationTrack_DocaMu2.push_back(-1);

            if(DocaMuon3.status()){ iIsolationTrack_DocaMu3.push_back(DocaMuon3.distance());}
            else iIsolationTrack_DocaMu3.push_back(-1);
            iIsolationTracksCollection.push_back(TheB->build(t));

            vector<TransientTrack> IsolationTrack_Muon1;
            vector<TransientTrack> IsolationTrack_Muon2;
            vector<TransientTrack> IsolationTrack_Muon3;


            IsolationTrack_Muon1.push_back(iTransientTracks.at(0));IsolationTrack_Muon1.push_back(TheB->build(t));
            IsolationTrack_Muon2.push_back(iTransientTracks.at(1));IsolationTrack_Muon2.push_back(TheB->build(t));
            IsolationTrack_Muon3.push_back(iTransientTracks.at(2));IsolationTrack_Muon3.push_back(TheB->build(t));

            KalmanVertexFitter KVF_IsoTrack_Mu1(true),KVF_IsoTrack_Mu2(true),KVF_IsoTrack_Mu3(true); 
            TransientVertex V_IsoTrack_Mu1,V_IsoTrack_Mu2,V_IsoTrack_Mu3;
            bool FitTrackMu1Ok(true);
            bool FitTrackMu2Ok(true);
            bool FitTrackMu3Ok(true);

            try {
               V_IsoTrack_Mu1=KVF_IsoTrack_Mu1.vertex(IsolationTrack_Muon1);
            } catch (...) {
               FitTrackMu1Ok = false;
            }

            try {
               V_IsoTrack_Mu2=KVF_IsoTrack_Mu2.vertex(IsolationTrack_Muon2);
            } catch (...) {
               FitTrackMu2Ok = false;
            }

            try {
               V_IsoTrack_Mu3=KVF_IsoTrack_Mu3.vertex(IsolationTrack_Muon3);
            } catch (...) {
               FitTrackMu3Ok = false;
            }

            iIsolationTrack_VertexWithSignalMuonIsValid.push_back(V_IsoTrack_Mu1.isValid());
            iIsolationTrack_VertexWithSignalMuonIsValid.push_back(V_IsoTrack_Mu2.isValid());
            iIsolationTrack_VertexWithSignalMuonIsValid.push_back(V_IsoTrack_Mu3.isValid());
            if(FitTrackMu1Ok&&FitTrackMu2Ok&&FitTrackMu3Ok);// dummy; to make it compile
            if(V_IsoTrack_Mu1.isValid())
            {
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(V_IsoTrack_Mu1.totalChiSquared());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu1.position().x());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu1.position().y());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu1.position().z());
            }
            else 
            {
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(-1);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
            }

            if(V_IsoTrack_Mu2.isValid())
            {	    
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(V_IsoTrack_Mu2.totalChiSquared());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu2.position().x());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu2.position().y());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu2.position().z());
            }
            else 
            {
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(-1);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
            }

            if(V_IsoTrack_Mu3.isValid())
            {	    
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(V_IsoTrack_Mu3.totalChiSquared());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu3.position().x());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu3.position().y());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu3.position().z());
            }
            else 
            {
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(-1);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
            }
            IsolationTrack_VertexWithSignalMuonIsValid.at(IsolationTrack_VertexWithSignalMuonIsValid.size() - 1).push_back(iIsolationTrack_VertexWithSignalMuonIsValid);
            IsolationTrack_VertexWithSignalMuonChi2.at(IsolationTrack_VertexWithSignalMuonChi2.size() - 1).push_back(iIsolationTrack_VertexWithSignalMuonChi2);
            IsolationTrack_VertexWithSignalMuonPosition.at(IsolationTrack_VertexWithSignalMuonPosition.size() - 1).push_back(iIsolationTrack_VertexWithSignalMuonPosition);
         }



         //--------------------------- Isolation Branch

         double dz = abs(t.dz(TheSecondaryVertexPoint));
         double dxy = abs(t.dxy(TheSecondaryVertexPoint));
         double dca_fv = sqrt(dz*dz+dxy*dxy);
         double dr_tau = deltaR(t.eta(), t.phi(), LVTau.Eta(), LVTau.Phi());

         ////////////////////////////////////////////////
         // Below is the isolation defined by Jian
         // iso no. 1b - using pt_min, drtau_max of the 3 mu
         if(t.pt() > 0.33*minmuon_pt && dr_tau < 3.*maxmuon_dr && dca_fv<0.05 ) {
            sumptalltracks += t.pt();
            sumalltracks++; // iso 3b
            if(dca_fv<mindist)mindist=dca_fv; // iso 4b
         } 

         if(t.pt()<1.0) continue;  // was 1.2
         // iso no. 1
         if(dr_tau < 0.5 && dca_fv<0.05 ) {
            sumptalltracks05 += t.pt();
            sumalltracks05++; // iso 3
            if(dca_fv<mindca_iso05)mindca_iso05=dca_fv; // iso 4
         }

         if(dca_fv<mindca_iso05)mindist05=dca_fv; // iso 4
         if(dca_fv<0.05)sumalltracks_b++; // iso 3b


         TransientTrack trkiso = theB->build(t);
         ClosestApproachInRPhi cAppm1, cAppm2, cAppm3;
         cAppm1.calculate(trkiso.initialFreeState(), iTransientTracks.at(0).initialFreeState());
         cAppm2.calculate(trkiso.initialFreeState(), iTransientTracks.at(1).initialFreeState());
         cAppm3.calculate(trkiso.initialFreeState(), iTransientTracks.at(2).initialFreeState());
         if(!(cAppm1.status()&&cAppm2.status()&&cAppm3.status())) continue;


         // iso no. 2
         if(deltaR(t.eta(), t.phi(), LV1.Eta(), LV1.Phi()) < 0.3 && cAppm1.distance() < 0.1) {// && dz1 < .3) 
            N_trk_1++;
            pt_trk_1 += t.pt();
         }
         if(deltaR(t.eta(), t.phi(), LV2.Eta(), LV2.Phi()) < 0.3 && cAppm2.distance() < 0.1) {//&& dz2 < .3) 
            N_trk_2++;
            pt_trk_2 += t.pt();
         }
         if(deltaR(t.eta(), t.phi(), LV3.Eta(), LV3.Phi()) < 0.3 && cAppm3.distance() < 0.1) {//&& dz3 < .3) 
            N_trk_3++;
            pt_trk_3 += t.pt();
         }
         if( (deltaR(t.eta(), t.phi(), LV1.Eta(), LV1.Phi()) < 0.3 && cAppm1.distance() < 0.1 )
               ||(deltaR(t.eta(), t.phi(), LV2.Eta(), LV2.Phi()) < 0.3 && cAppm2.distance() < 0.1 )
               ||(deltaR(t.eta(), t.phi(), LV3.Eta(), LV3.Phi()) < 0.3 && cAppm3.distance() < 0.1 )
           ) N_trk_total++;



         double dz_primaryvertex=abs(t.dz(pvPoint));

         if(!(dz_primaryvertex < 1))continue;
         double dxy_primaryvertex = abs(t.dxy(pvPoint));
         if(dxy_primaryvertex>0.1) N_trk0p1++;
         if(dxy_primaryvertex>0.2) N_trk0p2++;
         if(dxy_primaryvertex>0.5) N_trk0p5++;
         if(dxy_primaryvertex>maxdxy) maxdxy = dxy_primaryvertex;
      }
      //********************************************************************************
      //Here reconstruct the vertices of all signal candidates with all isolation tracks 


      IsolationTrack_dxySV.push_back(iIsolationTrack_dxySV);
      IsolationTrack_dzSV.push_back(iIsolationTrack_dzSV);

      IsolationTrack_dxyPV.push_back(iIsolationTrack_dxyPV);
      IsolationTrack_dzPV.push_back(iIsolationTrack_dzPV);
      IsolationTrack_DocaMu1.push_back(iIsolationTrack_DocaMu1);
      IsolationTrack_DocaMu2.push_back(iIsolationTrack_DocaMu2);
      IsolationTrack_DocaMu3.push_back(iIsolationTrack_DocaMu3);
      IsolationTrack_charge.push_back(iIsolationTrack_charge);
      IsolationTrack_isHighPurity.push_back(iIsolationTrack_isHighPurity);

      relative_iso = sumptalltracks/LVTau.Pt();
      relative_iso05 = sumptalltracks05/LVTau.Pt();
      relative_mu1_iso = pt_trk_1/LV1.Pt(); relative_mu2_iso = pt_trk_2/LV2.Pt(); relative_mu3_iso = pt_trk_3/LV3.Pt();
      relative_maxiso = TMath::Max(relative_mu1_iso, TMath::Max(relative_mu2_iso,relative_mu3_iso ));

      std::vector<float> isolation1, isolation2, isolation3, isolation4;
      isolation1.push_back(relative_iso);
      isolation1.push_back(sumalltracks);
      isolation1.push_back(mindist);

      isolation2.push_back(relative_iso05);
      isolation2.push_back(sumalltracks05);
      isolation2.push_back(mindist);
      isolation2.push_back(mindist05);


      isolation3.push_back(relative_mu1_iso);
      isolation3.push_back(relative_mu2_iso);
      isolation3.push_back(relative_mu3_iso);
      isolation3.push_back(relative_maxiso);

      isolation4.push_back(N_trk_1);
      isolation4.push_back(N_trk_2);
      isolation4.push_back(N_trk_3);

      isolation4.push_back(N_trk0p1);
      isolation4.push_back(N_trk0p2);
      isolation4.push_back(N_trk0p5);
      isolation4.push_back(maxdxy);

      Vertex_Isolation1.push_back(isolation1);
      Vertex_Isolation2.push_back(isolation2);
      Vertex_Isolation3.push_back(isolation3);
      Vertex_Isolation4.push_back(isolation4);
      index++;
   }

   for(size_t isv = 0; isv < svertices->size(); isv++) {
      const Vertex & sv = (*svertices)[isv];
      SV_Track_P4.push_back(std::vector<std::vector<float> >());
      std::vector<float> iSV_pos;
      iSV_pos.push_back(sv.x());
      iSV_pos.push_back(sv.y());
      iSV_pos.push_back(sv.z());
      SV_pos.push_back(iSV_pos);
      SV_Mass.push_back(sv.p4().M());

      TMatrixTSym<float> sv_cov(LorentzVectorParticle::NVertex);

      math::Error<3>::type sv_Cov;
      sv.fill(sv_Cov);

      for (int i = 0; i <LorentzVectorParticle::NVertex; i++){
         for (int j = 0; j < LorentzVectorParticle::NVertex; j++) {
            sv_cov(i, j) = sv_Cov(i, j);
            sv_cov(j, i) = sv_Cov(i, j);
         }
      }

      std::vector<float>  sv_covariance;
      for (int i = 0; i < LorentzVectorParticle::NVertex; i++) {
         for (int j = i; j < LorentzVectorParticle::NVertex; j++) {
            sv_covariance.push_back(sv_cov(i, j));
         }
      }
      SV_PosCovariance.push_back(sv_covariance);
      std::vector<int>    iSV_Trackcharge;
      for(Vertex::trackRef_iterator itk = sv.tracks_begin(); itk != sv.tracks_end(); itk++) {
         std::vector<float>  iSV_TrackP4;

         iSV_TrackP4.push_back(sqrt(pow((**itk).p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
         iSV_TrackP4.push_back((**itk).px());
         iSV_TrackP4.push_back((**itk).py());
         iSV_TrackP4.push_back((**itk).pz());
         SV_Track_P4.at(SV_Track_P4.size()-1).push_back(iSV_TrackP4);
         iSV_Trackcharge.push_back((**itk).charge());

      }
      SV_TrackCharge.push_back(iSV_Trackcharge);
   }
   return;
}

void T3MNtuple::fillVertices(const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<vector<pat::Muon> >& muons,
      const Handle<TrackCollection>& trackCollection, 
      const Handle<VertexCollection >& pvs,
      const Handle<std::vector<reco::VertexCompositePtrCandidate> >& svertices,
      const Handle<reco::BeamSpot>& beamSpotHandle)
{

   const BeamSpot& bs = *beamSpotHandle;

   std::vector<std::vector<std::vector<double> > > particles_p4;
   std::vector<std::vector<TransientTrack> > signalTracksCollection;
   ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
   if(ThreeMuons_idx.size()!=0){
      for ( auto &iThreeMuon :  ThreeMuons_idx ) {
         particles_p4.push_back(std::vector<std::vector<double> > ());
         vector<TransientTrack> isignalTracksCollection;
         ESHandle<TransientTrackBuilder> theB;
         iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
         for ( auto &iMuon :  iThreeMuon ) {
            pat::MuonRef Muon(muons, iMuon);

            //edm::View<reco::Muon>::const_reference Muon = muons->at(iMuon);

            TrackRef MuonTrack = Muon->innerTrack();
            isignalTracksCollection.push_back(theB->build(MuonTrack));
            std::vector<double> iiparticles_p4;
            iiparticles_p4.push_back(Muon->p4().E());
            iiparticles_p4.push_back(Muon->p4().Px());
            iiparticles_p4.push_back(Muon->p4().Py());
            iiparticles_p4.push_back(Muon->p4().Pz());
            particles_p4.at(particles_p4.size() -1).push_back(iiparticles_p4);
         }
         signalTracksCollection.push_back(isignalTracksCollection);
      }
   }

   if(TwoMuonsTrack_idx.size()!=0){
      std::vector<std::vector<double> > iparticle_p4;
      for ( auto &iTwoMuonsTracks :  TwoMuonsTrack_idx ) {

         vector<TransientTrack> isignalTracksCollection;

         pat::MuonRef Muon1(muons, iTwoMuonsTracks.at(0));
         pat::MuonRef Muon2(muons, iTwoMuonsTracks.at(1));

         TrackRef track1 = Muon1->innerTrack();
         TrackRef track2 = Muon2->innerTrack();
         TrackRef track3 = TrackRef(trackCollection, iTwoMuonsTracks.at(2));

         isignalTracksCollection.push_back(theB->build(track1));
         isignalTracksCollection.push_back(theB->build(track2));
         isignalTracksCollection.push_back(theB->build(track3));
         signalTracksCollection.push_back(isignalTracksCollection);

         std::vector<double> particle1_p4;
         std::vector<double> particle2_p4;
         std::vector<double> particle3_p4;

         particle1_p4.push_back(Muon1->p4().E());       particle2_p4.push_back(Muon2->p4().E());       particle3_p4.push_back(sqrt(pow(track3->p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
         particle1_p4.push_back(Muon1->p4().Px());      particle2_p4.push_back(Muon2->p4().Px());      particle3_p4.push_back(track3->px());
         particle1_p4.push_back(Muon1->p4().Py());      particle2_p4.push_back(Muon2->p4().Py());      particle3_p4.push_back(track3->py());
         particle1_p4.push_back(Muon1->p4().Pz());      particle2_p4.push_back(Muon2->p4().Pz());      particle3_p4.push_back(track3->pz());
         iparticle_p4.push_back(particle1_p4);          iparticle_p4.push_back(particle2_p4);          iparticle_p4.push_back(particle3_p4);
         particles_p4.push_back(iparticle_p4);
      }
   }

   if (DEBUG)   std::cout<<"size  "<< signalTracksCollection.size() << std::endl;
   unsigned int index(0);
   for ( auto &iTransientTracks :  signalTracksCollection ){
      Vertex_signal_KF_pos.push_back(std::vector<double> ());
      Vertex_signal_KF_cov.push_back(std::vector<double> ());
      Vertex_signal_KF_refittedTracksP4.push_back(std::vector<std::vector<double> >());

      Vertex_signal_AF_pos.push_back(std::vector<double> ());
      ClosestApproachInRPhi cApp12, cApp23, cApp31;
      cApp12.calculate(iTransientTracks[0].initialFreeState(), iTransientTracks[1].initialFreeState());
      cApp23.calculate(iTransientTracks[1].initialFreeState(), iTransientTracks[2].initialFreeState());
      cApp31.calculate(iTransientTracks[2].initialFreeState(), iTransientTracks[0].initialFreeState());
      std::vector<double> iVertex_signal_dca_reco;
      if((cApp12.status()&&cApp23.status()&&cApp31.status())) { 
         iVertex_signal_dca_reco.push_back(cApp12.distance());   // order 12,23,31, max
         iVertex_signal_dca_reco.push_back(cApp23.distance());
         iVertex_signal_dca_reco.push_back(cApp31.distance());
         iVertex_signal_dca_reco.push_back(TMath::Max(dca12_reco, TMath::Max(dca23_reco, dca31_reco)));
      } else {
         iVertex_signal_dca_reco.push_back(-1.);
         iVertex_signal_dca_reco.push_back(-1.);
         iVertex_signal_dca_reco.push_back(-1.);
         iVertex_signal_dca_reco.push_back(-1.);
      }
      Vertex_signal_dca_reco.push_back(iVertex_signal_dca_reco);

      TransientVertex transVtx;
      KalmanVertexFitter kvf(true);
      bool FitOk(true);
      try {
         transVtx = kvf.vertex(iTransientTracks); 
      } catch (...) {
         FitOk = false;
      }
      if (!transVtx.hasRefittedTracks())
         FitOk = false;
      if (transVtx.refittedTracks().size() != iTransientTracks.size())
         FitOk = false;
      TLorentzVector ThreeCandidate(0,0,0,0);

      math::XYZPoint TheSecondaryVertexPoint;
      if(FitOk){
         Vertex_signal_KF_Chi2.push_back(transVtx.totalChiSquared());
         Vertex_signal_KF_pos.at(index).push_back(transVtx.position().x());
         Vertex_signal_KF_pos.at(index).push_back(transVtx.position().y());
         Vertex_signal_KF_pos.at(index).push_back(transVtx.position().z());

         reco::Vertex secondaryVertex = transVtx;
         TMatrixTSym<double> svcov(LorentzVectorParticle::NVertex);
         math::Error<3>::type svCov;
         secondaryVertex.fill(svCov);

         for (int i = 0; i <LorentzVectorParticle::NVertex; i++){
            for (int j = 0; j < LorentzVectorParticle::NVertex; j++) {
               svcov(i, j) = svCov(i, j);
               svcov(j, i) = svCov(i, j);
            }
         }
         for (int i = 0; i < LorentzVectorParticle::NVertex; i++) {
            for (int j = i; j < LorentzVectorParticle::NVertex; j++) {
               Vertex_signal_KF_cov.at(index).push_back(svcov(i, j));
            }
         }

         TheSecondaryVertexPoint = math::XYZPoint(transVtx.position().x(), transVtx.position().y(), transVtx.position().z());
         vector<TransientTrack>::const_iterator trkIt = transVtx.refittedTracks().begin();

         reco::Vertex TripletVtx = reco::Vertex(TheSecondaryVertexPoint, svCov, transVtx.totalChiSquared(), 3, 3);

         VertexState BSstate(bs);
         VertexDistanceXY vertTool2D;
         double BSdistance2D = vertTool2D.distance(BSstate, TripletVtx).value();
         double BSdist_err2D = vertTool2D.distance(BSstate, TripletVtx).error();
         double BSdist_sign2D =vertTool2D.distance(BSstate, TripletVtx).significance();

         Vertex_signal_KF_BS_2Ddistance.push_back(BSdistance2D);
         Vertex_signal_KF_BS_error.push_back(BSdist_err2D);
         Vertex_signal_KF_BS_significance.push_back(BSdist_sign2D);

         for(; trkIt != transVtx.refittedTracks().end(); ++ trkIt) {
            std::vector<double> irefitted_tracks_p4;
            const Track & trkrefit = trkIt->track();
            irefitted_tracks_p4.push_back(sqrt(pow(trkrefit.p(),2.0) + pow(PDGInfo::mu_mass(),2.0)));
            irefitted_tracks_p4.push_back(trkrefit.px());
            irefitted_tracks_p4.push_back(trkrefit.py());
            irefitted_tracks_p4.push_back(trkrefit.pz());
            ThreeCandidate+=TLorentzVector(trkrefit.px(),trkrefit.py(),trkrefit.pz(),sqrt(pow(trkrefit.p(),2.0) + pow(PDGInfo::mu_mass(),2.0)));
            Vertex_signal_KF_refittedTracksP4.at(Vertex_signal_KF_refittedTracksP4.size() -1).push_back(irefitted_tracks_p4);
         }
      }

      bool AFitOk(true);
      AdaptiveVertexFitter avf;
      TransientVertex AdaptivetransVtx;
      avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
      try {
         AdaptivetransVtx = avf.vertex(iTransientTracks);
      } catch (...) {
         AFitOk = false;
      }
      if(!AdaptivetransVtx.isValid()) 
         AFitOk = false;
      if (AFitOk){ 
         Vertex_signal_AF_Chi2.push_back(AdaptivetransVtx.totalChiSquared());
         Vertex_signal_AF_Ndf.push_back(AdaptivetransVtx.degreesOfFreedom());
         Vertex_signal_AF_pos.at(index).push_back(AdaptivetransVtx.position().x());
         Vertex_signal_AF_pos.at(index).push_back(AdaptivetransVtx.position().y());
         Vertex_signal_AF_pos.at(index).push_back(AdaptivetransVtx.position().z());

      } else {
         Vertex_signal_AF_Chi2.push_back(-1);
      }

      //--------------------  Fit each track pair 
      vector<TransientTrack> trackpair12, trackpair23, trackpair31;
      trackpair12.push_back(iTransientTracks.at(0)); trackpair12.push_back(iTransientTracks.at(1));
      trackpair23.push_back(iTransientTracks.at(1)); trackpair23.push_back(iTransientTracks.at(2));
      trackpair31.push_back(iTransientTracks.at(2)); trackpair31.push_back(iTransientTracks.at(0));
      KalmanVertexFitter kvf_trks12(true), kvf_trks23(true), kvf_trks31(true);
      TransientVertex fv_trks12;// = kvf_trks12.vertex(trackpair12);
      TransientVertex fv_trks23;// = kvf_trks23.vertex(trackpair23);
      TransientVertex fv_trks31;// = kvf_trks31.vertex(trackpair31);


      bool Fit1Ok(true);
      try {
         fv_trks12 = kvf_trks12.vertex(trackpair12); 
      } catch (...) {
         Fit1Ok = false;
      }

      bool Fit2Ok(true);
      try {
         fv_trks23 = kvf_trks23.vertex(trackpair23); 
      } catch (...) {
         Fit2Ok = false;
      }

      bool Fit3Ok(true);
      try {
         fv_trks31 = kvf_trks31.vertex(trackpair31); 
      } catch (...) {
         Fit3Ok = false;
      }

      std::vector<double> iVertex_pair_quality;
      std::vector<double> iVertex_pairfit_status;
      iVertex_pairfit_status.push_back(Fit1Ok);
      iVertex_pairfit_status.push_back(Fit2Ok);
      iVertex_pairfit_status.push_back(Fit3Ok);

      std::vector<float> iVertex_Pair12_Pos;
      std::vector<float> iVertex_Pair23_Pos;
      std::vector<float> iVertex_Pair31_Pos;
      if(fv_trks12.isValid())
      {
         iVertex_pair_quality.push_back(fv_trks12.totalChiSquared());
         iVertex_Pair12_Pos.push_back( fv_trks12.position().x());
         iVertex_Pair12_Pos.push_back( fv_trks12.position().y());
         iVertex_Pair12_Pos.push_back( fv_trks12.position().z());
      }

      else
      {
         iVertex_pair_quality.push_back(-1);
         iVertex_Pair12_Pos.push_back( 99.);
         iVertex_Pair12_Pos.push_back( 99.);
         iVertex_Pair12_Pos.push_back( 99.);
      }
      if(fv_trks23.isValid())
      {
         iVertex_pair_quality.push_back(fv_trks23.totalChiSquared());
         iVertex_Pair23_Pos.push_back( fv_trks23.position().x());
         iVertex_Pair23_Pos.push_back( fv_trks23.position().y());
         iVertex_Pair23_Pos.push_back( fv_trks23.position().z());
      }
      else
      {
         iVertex_pair_quality.push_back(-1);
         iVertex_Pair23_Pos.push_back( 99.);
         iVertex_Pair23_Pos.push_back( 99.);
         iVertex_Pair23_Pos.push_back( 99.);
      }
      if(fv_trks31.isValid())
      {
         iVertex_pair_quality.push_back(fv_trks31.totalChiSquared());
         iVertex_Pair31_Pos.push_back( fv_trks31.position().x());
         iVertex_Pair31_Pos.push_back( fv_trks31.position().y());
         iVertex_Pair31_Pos.push_back( fv_trks31.position().z());
      }
      else
      {
         iVertex_pair_quality.push_back(-1);
         iVertex_Pair31_Pos.push_back( 99.);
         iVertex_Pair31_Pos.push_back( 99.);
         iVertex_Pair31_Pos.push_back( 99.);
      }

      Vertex_Pair12_Pos.push_back( iVertex_Pair12_Pos);
      Vertex_Pair23_Pos.push_back( iVertex_Pair23_Pos);
      Vertex_Pair31_Pos.push_back( iVertex_Pair31_Pos);


      Vertex_pair_quality.push_back(iVertex_pair_quality);
      Vertex_pairfit_status.push_back(iVertex_pairfit_status);

      ///////////////////////////////////////////
      //  find here the primary vertex with the best
      //  alignement to the tri-muon 
      double dphi_pv = -1.0;
      unsigned int primaryvertex_index;
      for(unsigned int vertex_index = 0; vertex_index  < pvs->size(); vertex_index++) {
         const Vertex & pvertex = (*pvs)[vertex_index];
         TVector3 Dv3D_reco(transVtx.position().x() - pvertex.x(), transVtx.position().y() - pvertex.y(), transVtx.position().z() - pvertex.z());
         double Cosdphi_3D = Dv3D_reco.Dot(ThreeCandidate.Vect())/(Dv3D_reco.Mag()*ThreeCandidate.Vect().Mag());
         if(Cosdphi_3D>dphi_pv){
            dphi_pv = Cosdphi_3D;
            primaryvertex_index=vertex_index;
         }
      }

      const Vertex & MatchedPrimaryVertex = (*pvs)[primaryvertex_index];
      dump_pv_index_to_fill.push_back(primaryvertex_index);


      std::vector<double>  iprimaryVertex_Pos;
      iprimaryVertex_Pos.push_back(MatchedPrimaryVertex.x());
      iprimaryVertex_Pos.push_back(MatchedPrimaryVertex.y());
      iprimaryVertex_Pos.push_back(MatchedPrimaryVertex.z());
      Vertex_MatchedPrimaryVertex.push_back(iprimaryVertex_Pos);

      double tempdz(99.);
      unsigned int secondbest_primaryvertex_index(0);
      for(unsigned int vertex_index = 0; vertex_index  < pvs->size(); vertex_index++) {
         if(vertex_index == primaryvertex_index) continue;
         const Vertex & temp_pv = (*pvs)[vertex_index];
         if(fabs(temp_pv.z() -  MatchedPrimaryVertex.z()) < tempdz ){
            tempdz = fabs(temp_pv.z() -  MatchedPrimaryVertex.z());
            secondbest_primaryvertex_index = vertex_index;
         }
      }


      const Vertex & SecondBestPrimaryVertex = (*pvs)[secondbest_primaryvertex_index];
      std::vector<double> iSecondBestprimaryVertex_Pos;
      if( pvs->size()>1){
         iSecondBestprimaryVertex_Pos.push_back(SecondBestPrimaryVertex.x());
         iSecondBestprimaryVertex_Pos.push_back(SecondBestPrimaryVertex.y());
         iSecondBestprimaryVertex_Pos.push_back(SecondBestPrimaryVertex.z());
      }else{
         iSecondBestprimaryVertex_Pos.push_back(-99.);
         iSecondBestprimaryVertex_Pos.push_back(-99.);
         iSecondBestprimaryVertex_Pos.push_back(-99.);
      }
      Vertex_SecondBestPrimaryVertex.push_back(iSecondBestprimaryVertex_Pos);


      int NParticlesComingFromPV(0);
      vector<TransientTrack> primaryvertexTransientTracks;// remove muon candidates from the PV to perform refit
      for(Vertex::trackRef_iterator itk = MatchedPrimaryVertex.tracks_begin(); itk != MatchedPrimaryVertex.tracks_end(); itk++) {
         if((**itk).pt()>1) {
            if(deltaR(iTransientTracks.at(0).track().eta(), iTransientTracks.at(0).track().phi(), (**itk).eta(), (**itk).phi())<0.01){NParticlesComingFromPV++;continue;}
            if(deltaR(iTransientTracks.at(1).track().eta(), iTransientTracks.at(1).track().phi(), (**itk).eta(), (**itk).phi())<0.01){NParticlesComingFromPV++;continue;}
            if(deltaR(iTransientTracks.at(2).track().eta(), iTransientTracks.at(2).track().phi(), (**itk).eta(), (**itk).phi())<0.01){NParticlesComingFromPV++;continue;}
         }
         primaryvertexTransientTracks.push_back(theB->build(**itk));
      }
      Vertex_NMuonsAssocWithPV.push_back(NParticlesComingFromPV);
      KalmanVertexFitter pv_fit(true);
      bool FitPVOk(true);
      TransientVertex pvvertex;

      if(primaryvertexTransientTracks.size() >1){
         try {
            pvvertex =  pv_fit.vertex(primaryvertexTransientTracks);
         } catch (...) {
            FitPVOk = false;
         }
      }

      Vertex_RefitPVisValid.push_back(pvvertex.isValid());
      std::vector<double> iRefitprimaryVertex_Pos;
      if(FitPVOk && pvvertex.isValid()){
         iRefitprimaryVertex_Pos.push_back(pvvertex.position().x());
         iRefitprimaryVertex_Pos.push_back(pvvertex.position().y());
         iRefitprimaryVertex_Pos.push_back(pvvertex.position().z());
      }
      Vertex_MatchedRefitPrimaryVertex.push_back(iRefitprimaryVertex_Pos);


      TMatrixTSym<double> pvcov(LorentzVectorParticle::NVertex);

      math::Error<3>::type pvCov;
      MatchedPrimaryVertex.fill(pvCov);

      for (int i = 0; i <LorentzVectorParticle::NVertex; i++){
         for (int j = 0; j < LorentzVectorParticle::NVertex; j++) {
            pvcov(i, j) = pvCov(i, j);
            pvcov(j, i) = pvCov(i, j);
         }
      }

      std::vector<double>  pv_cov;     
      for (int i = 0; i < LorentzVectorParticle::NVertex; i++) {
         for (int j = i; j < LorentzVectorParticle::NVertex; j++) {
            pv_cov.push_back(pvcov(i, j));
         }
      }

      Vertex_MatchedRefitPrimaryVertex_covariance.push_back(pv_cov);

      Vertex final_pv = MatchedPrimaryVertex;  
      if(pvvertex.isValid()) final_pv = Vertex(pvvertex);

      math::XYZPoint pvPoint = math::XYZPoint(final_pv.x(), final_pv.y(), final_pv.z());
      math::XYZPoint bsPoint = math::XYZPoint(beamSpotHandle->position().x(), beamSpotHandle->position().y(), beamSpotHandle->position().z());



      std::vector<double> iVertex_d0BeamSpot_reco;
      iVertex_d0BeamSpot_reco.push_back(abs(iTransientTracks.at(0).track().dxy(bsPoint)));
      iVertex_d0BeamSpot_reco.push_back(abs(iTransientTracks.at(1).track().dxy(bsPoint)));
      iVertex_d0BeamSpot_reco.push_back(abs(iTransientTracks.at(2).track().dxy(bsPoint)));
      Vertex_d0BeamSpot_reco.push_back(iVertex_d0BeamSpot_reco);


      std::vector<double> iVertex_d0BeamSpot_reco_sig;
      double d0ErrorToBs_1  = sqrt( iTransientTracks.at(0).track().d0Error() * iTransientTracks.at(0).track().d0Error() +
            0.5*  beamSpotHandle->BeamWidthX()* beamSpotHandle->BeamWidthX()+
            0.5*  beamSpotHandle->BeamWidthY()* beamSpotHandle->BeamWidthY() );

      double d0ErrorToBs_2  = sqrt( iTransientTracks.at(1).track().d0Error() * iTransientTracks.at(1).track().d0Error() +
            0.5*  beamSpotHandle->BeamWidthX()* beamSpotHandle->BeamWidthX()+
            0.5*  beamSpotHandle->BeamWidthY()* beamSpotHandle->BeamWidthY() );

      double d0ErrorToBs_3  = sqrt( iTransientTracks.at(2).track().d0Error() * iTransientTracks.at(2).track().d0Error() +
            0.5*  beamSpotHandle->BeamWidthX()* beamSpotHandle->BeamWidthX()+
            0.5*  beamSpotHandle->BeamWidthY()* beamSpotHandle->BeamWidthY() );



      if(d0ErrorToBs_1!=0){  iVertex_d0BeamSpot_reco_sig.push_back( abs(iTransientTracks.at(0).track().dxy(bsPoint)) / d0ErrorToBs_1);} else {iVertex_d0BeamSpot_reco_sig.push_back(-1);}
      if(d0ErrorToBs_2!=0){  iVertex_d0BeamSpot_reco_sig.push_back( abs(iTransientTracks.at(1).track().dxy(bsPoint)) / d0ErrorToBs_2);} else {iVertex_d0BeamSpot_reco_sig.push_back(-1);}
      if(d0ErrorToBs_3!=0){  iVertex_d0BeamSpot_reco_sig.push_back( abs(iTransientTracks.at(1).track().dxy(bsPoint)) / d0ErrorToBs_3);} else {iVertex_d0BeamSpot_reco_sig.push_back(-1);}

      Vertex_d0BeamSpot_reco_sig.push_back(iVertex_d0BeamSpot_reco_sig);


      std::vector<double> iVertex_d0SV_reco;
      iVertex_d0SV_reco.push_back(abs(iTransientTracks.at(0).track().dxy(TheSecondaryVertexPoint)));
      iVertex_d0SV_reco.push_back(abs(iTransientTracks.at(1).track().dxy(TheSecondaryVertexPoint)));
      iVertex_d0SV_reco.push_back(abs(iTransientTracks.at(2).track().dxy(TheSecondaryVertexPoint)));
      Vertex_d0SV_reco.push_back(iVertex_d0SV_reco);


      std::vector<double> iVertex_dzSV_reco;
      iVertex_dzSV_reco.push_back(abs(iTransientTracks.at(0).track().dz(TheSecondaryVertexPoint)));
      iVertex_dzSV_reco.push_back(abs(iTransientTracks.at(1).track().dz(TheSecondaryVertexPoint)));
      iVertex_dzSV_reco.push_back(abs(iTransientTracks.at(2).track().dz(TheSecondaryVertexPoint)));
      Vertex_dzSV_reco.push_back(iVertex_dzSV_reco);



      std::vector<double> iVertex_d0_reco;
      iVertex_d0_reco.push_back(abs(iTransientTracks.at(0).track().dxy(pvPoint)));
      iVertex_d0_reco.push_back(abs(iTransientTracks.at(1).track().dxy(pvPoint)));
      iVertex_d0_reco.push_back(abs(iTransientTracks.at(2).track().dxy(pvPoint)));
      Vertex_d0_reco.push_back(iVertex_d0_reco);


      std::vector<double> iVertex_dz_reco;
      iVertex_dz_reco.push_back(abs(iTransientTracks.at(0).track().dz(pvPoint)));
      iVertex_dz_reco.push_back(abs(iTransientTracks.at(1).track().dz(pvPoint)));
      iVertex_dz_reco.push_back(abs(iTransientTracks.at(2).track().dz(pvPoint)));
      Vertex_dz_reco.push_back(iVertex_dz_reco);



      TLorentzVector LV1=TLorentzVector(particles_p4.at(index).at(0).at(1), particles_p4.at(index).at(0).at(2), particles_p4.at(index).at(0).at(3),particles_p4.at(index).at(0).at(0));
      TLorentzVector LV2=TLorentzVector(particles_p4.at(index).at(1).at(1), particles_p4.at(index).at(1).at(2), particles_p4.at(index).at(1).at(3),particles_p4.at(index).at(1).at(0));
      TLorentzVector LV3=TLorentzVector(particles_p4.at(index).at(2).at(1), particles_p4.at(index).at(2).at(2), particles_p4.at(index).at(2).at(3),particles_p4.at(index).at(2).at(0));
      TLorentzVector LVTau = LV1 + LV2 + LV3;
      m3mu_reco = LVTau.M();

      GlobalVector dir1(particles_p4.at(index).at(0).at(1), particles_p4.at(index).at(0).at(2), particles_p4.at(index).at(0).at(3));
      GlobalVector dir2(particles_p4.at(index).at(1).at(1), particles_p4.at(index).at(1).at(2), particles_p4.at(index).at(1).at(3));
      GlobalVector dir3(particles_p4.at(index).at(2).at(1), particles_p4.at(index).at(2).at(2), particles_p4.at(index).at(2).at(3));

      std::pair<bool, Measurement1D> ip2d_1 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(0), dir1, final_pv);
      std::pair<bool, Measurement1D> ip2d_2 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(1), dir2, final_pv);
      std::pair<bool, Measurement1D> ip2d_3 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(2), dir3, final_pv);


      std::pair<bool, Measurement1D> ip2dSV_1 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(0), dir1, transVtx);
      std::pair<bool, Measurement1D> ip2dSV_2 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(1), dir2, transVtx);
      std::pair<bool, Measurement1D> ip2dSV_3 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(2), dir3, transVtx);


      //  std::pair<bool, Measurement1D> ipBS2d_1 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(0), dir1, );
      //  std::pair<bool, Measurement1D> ipBS2d_2 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(1), dir2, final_pv);
      //  std::pair<bool, Measurement1D> ipBS2d_3 = IPTools::signedTransverseImpactParameter(iTransientTracks.at(2), dir3, final_pv);

      std::vector<double> iVertex_d0sig_reco;
      if(ip2d_1.first){ iVertex_d0sig_reco.push_back(abs(ip2d_1.second.value()/ip2d_1.second.error()));} else{iVertex_d0sig_reco.push_back(-1);}
      if(ip2d_2.first){ iVertex_d0sig_reco.push_back(abs(ip2d_2.second.value()/ip2d_2.second.error()));} else{iVertex_d0sig_reco.push_back(-1);}
      if(ip2d_3.first){ iVertex_d0sig_reco.push_back(abs(ip2d_3.second.value()/ip2d_3.second.error()));} else{iVertex_d0sig_reco.push_back(-1);}
      Vertex_d0sig_reco.push_back(iVertex_d0sig_reco);


      std::vector<double> iVertex_d0sigSV_reco;
      if(ip2dSV_1.first){ iVertex_d0sigSV_reco.push_back(abs(ip2dSV_1.second.value()/ip2dSV_1.second.error()));} else{iVertex_d0sigSV_reco.push_back(-1);}
      if(ip2dSV_2.first){ iVertex_d0sigSV_reco.push_back(abs(ip2dSV_2.second.value()/ip2dSV_2.second.error()));} else{iVertex_d0sigSV_reco.push_back(-1);}
      if(ip2dSV_3.first){ iVertex_d0sigSV_reco.push_back(abs(ip2dSV_3.second.value()/ip2dSV_3.second.error()));} else{iVertex_d0sigSV_reco.push_back(-1);}
      Vertex_d0sigSV_reco.push_back(iVertex_d0sigSV_reco);


      TVector3 dv2D_reco(-final_pv.position().x() + TheSecondaryVertexPoint.x(), -final_pv.position().y() + TheSecondaryVertexPoint.y(), 0);
      TVector3 vtauxy(ThreeCandidate.Px(), ThreeCandidate.Py(), 0);
      fv_cosdphi = dv2D_reco.Dot(vtauxy)/(dv2D_reco.Perp()*vtauxy.Perp());
      VertexDistanceXY vdistXY;
      Measurement1D distXY = vdistXY.distance(Vertex(transVtx), final_pv);
      std::vector<double> iVertex_2Ddisplacement;
      iVertex_2Ddisplacement.push_back(distXY.value());
      iVertex_2Ddisplacement.push_back(distXY.significance());
      iVertex_2Ddisplacement.push_back(distXY.value()*fv_cosdphi * m3mu_reco/vtauxy.Perp());
      Vertex_2Ddisplacement.push_back(iVertex_2Ddisplacement);


      TVector3 vtauxyz(ThreeCandidate.Px(), ThreeCandidate.Py(), ThreeCandidate.Pz());
      TVector3 dv3D_reco(-final_pv.position().x() + TheSecondaryVertexPoint.x(), -final_pv.position().y() + TheSecondaryVertexPoint.y(), -final_pv.position().z() + TheSecondaryVertexPoint.z());
      fv_cosdphi3D = dv3D_reco.Dot(vtauxyz)/(dv3D_reco.Mag()*vtauxyz.Mag());
      VertexDistance3D dist;

      std::vector<double> iVertex_3Ddisplacement;
      iVertex_3Ddisplacement.push_back(dist.distance(Vertex(transVtx), final_pv).value());
      iVertex_3Ddisplacement.push_back(dist.distance(Vertex(transVtx), final_pv).significance());
      iVertex_3Ddisplacement.push_back(fv_d3D*fv_cosdphi3D*m3mu_reco/ThreeCandidate.P());
      Vertex_3Ddisplacement.push_back(iVertex_3Ddisplacement);


      ///////////////////////////////////////////////////////////////
      //    Here fill the isolation


      IsolationBranch_Trackp4.push_back(std::vector<std::vector<float> >());
      if (DEBUG) cout << "Size of the nTracks matched to PV: "<<MatchedPrimaryVertex.tracksSize()<<endl;
      //for(Vertex::trackRef_iterator itk = MatchedPrimaryVertex.tracks_begin(); itk != MatchedPrimaryVertex.tracks_end(); itk++) {}
      for (unsigned int iTrack = 0; iTrack<(*trackCollection).size(); iTrack++){

         TrackRef track = reco::TrackRef(trackCollection, iTrack);
         TransientTrack trans_track = theB->build(track);
         
         // filter tracks (same definition as used in /RecoVertex/PrimaryVertexProducer/python/OfflinePrimaryVertices_cfi.py in CMSSW)
         bool IPSigCut = (trans_track.stateAtBeamLine().transverseImpactParameter().significance() < 4.0) &&
                         (trans_track.stateAtBeamLine().transverseImpactParameter().error() < 1.0) &&
                         (trans_track.track().dzError() < 1.0);
         if (!IPSigCut) continue;
         if (abs(trans_track.impactPointState().globalMomentum().eta())>2.4) continue;
         if (trans_track.normalizedChi2()>10.0) continue;
         if (trans_track.hitPattern().pixelLayersWithMeasurement()<2) continue;
         if (trans_track.hitPattern().trackerLayersWithMeasurement()<5) continue;
         if (trans_track.track().quality(reco::TrackBase::undefQuality)) continue; // returns loose track by default

         // Additional requirements ( diaplcement from PV)
         if (abs(trans_track.track().dz(MatchedPrimaryVertex.position()))>0.5) continue;
         if (abs(trans_track.track().dxy(MatchedPrimaryVertex.position()))>0.2) continue;

         if(deltaR(iTransientTracks.at(0).track().eta(), iTransientTracks.at(0).track().phi(), track->eta(), track->phi())<0.001)continue;
         if(deltaR(iTransientTracks.at(1).track().eta(), iTransientTracks.at(1).track().phi(), track->eta(), track->phi())<0.001)continue;
         if(deltaR(iTransientTracks.at(2).track().eta(), iTransientTracks.at(2).track().phi(), track->eta(), track->phi())<0.001)continue;

         std::vector<float> iIsolationBranch_Track_p4;

         iIsolationBranch_Track_p4.push_back(sqrt(pow(track->p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
         iIsolationBranch_Track_p4.push_back(track->px());
         iIsolationBranch_Track_p4.push_back(track->py());
         iIsolationBranch_Track_p4.push_back(track->pz());

         IsolationBranch_Trackp4.at(IsolationBranch_Trackp4.size() - 1).push_back(iIsolationBranch_Track_p4);
      }
         if (DEBUG) cout<<"Size of nTracks found close to PV: "<<IsolationBranch_Trackp4.at(IsolationBranch_Trackp4.size()-1).size()<<endl;



      //  former Isolation
      float minmuon_pt(999.), maxmuon_dr(0.);

      if(LV1.Pt()<minmuon_pt)minmuon_pt=LV1.Pt();
      if(LV2.Pt()<minmuon_pt)minmuon_pt=LV2.Pt();
      if(LV3.Pt()<minmuon_pt)minmuon_pt=LV3.Pt();

      float  drLV1Tau(deltaR(LV1.Eta(), LV1.Phi(), LVTau.Eta(), LVTau.Phi()));
      float  drLV2Tau(deltaR(LV2.Eta(), LV2.Phi(), LVTau.Eta(), LVTau.Phi()));
      float  drLV3Tau(deltaR(LV3.Eta(), LV3.Phi(), LVTau.Eta(), LVTau.Phi()));

      maxmuon_dr = TMath::Max(drLV1Tau, TMath::Max(drLV2Tau,drLV3Tau));

      float sumptalltracks(0.), sumalltracks(0.), mindist(99.);
      float sumptalltracks05(0.), sumalltracks05(0.), mindist05(99.), sumalltracks_b(0.);
      float pt_trk_1(0.), pt_trk_2(0.), pt_trk_3(0.);
      float N_trk_1(0.),N_trk_2(0.),N_trk_3(0.), N_trk_total(0.);
      float N_trk0p1(0.), N_trk0p2(0.), N_trk0p5(0.), maxdxy(0.);
      float relative_iso(0.), relative_iso05(0.), relative_mu1_iso(0.),relative_mu2_iso(0.),relative_mu3_iso(0.);
      float relative_maxiso(0.);


      //--------------------------- Isolation Branch

      std::vector<TLorentzVector>  MuLVs;
      MuLVs.push_back(LV1);
      MuLVs.push_back(LV2);
      MuLVs.push_back(LV3);


      std::vector<int> sortedindices = SortByPt(MuLVs);


      IsolationTrack_p4.push_back(std::vector<std::vector<float> >());
      IsolationTrack_VertexWithSignalMuonIsValid.push_back(std::vector<std::vector<int> >());
      IsolationTrack_VertexWithSignalMuonChi2.push_back(std::vector<std::vector<float> >());
      IsolationTrack_VertexWithSignalMuonPosition.push_back(std::vector<std::vector<float> >());

      std::vector<float> iIsolationTrack_dxySV;
      std::vector<float> iIsolationTrack_dzSV;
      std::vector<float> iIsolationTrack_dxyPV;
      std::vector<float> iIsolationTrack_dzPV;
      std::vector<float> iIsolationTrack_DocaMu1;
      std::vector<float> iIsolationTrack_DocaMu2;
      std::vector<float> iIsolationTrack_DocaMu3;
      std::vector<int>   iIsolationTrack_charge;
      std::vector<int>   iIsolationTrack_isHighPurity;
      std::vector<int>   IsolationTracksIndices;


      vector<TransientTrack> iIsolationTracksCollection;
      ESHandle<TransientTrackBuilder> TheB;
      std::vector<std::vector<int> >  VertexWithSignalMuonIsValid;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TheB);



      if (DEBUG)     std::cout<<" sig   " << ThreeMuons_idx.size() << " two mu plus track   "<< TwoMuonsTrack_idx.size() <<  " signaltrack collect     "<< signalTracksCollection.size() <<  " signal Counter  "<<index <<std::endl;

      for(size_t iTrack = 0; iTrack < (*trackCollection).size(); iTrack++) {

         TrackRef t = reco::TrackRef(trackCollection, iTrack);


         if(!(t->quality(TrackBase::highPurity)))continue; //-- this might be weaker
         if (DEBUG) cout<<"Found tight track"<<endl;
         if(deltaR(LV1.Eta(), LV1.Phi(), t->eta(), t->phi())<0.01)continue;
         if(deltaR(LV2.Eta(), LV2.Phi(), t->eta(), t->phi())<0.01)continue;
         if(deltaR(LV3.Eta(), LV3.Phi(), t->eta(), t->phi())<0.01)continue;
            if (DEBUG) cout<<"Track not matched to muons"<<endl;
         if(index>=ThreeMuons_idx.size()) break;  // <----- fillinf isolation vertexing only for a signal candidate
         //if(abs(t.dz(pvPoint))< 0.5 && t.quality(TrackBase::tight) && sqrt(t.px()*t.px() + t.py()*t.py() ) > 0.5){//  && deltaR(t.eta(), t.phi(), LVTau.Eta(), LVTau.Phi()) < 1.){ }}
         if(abs(t->dz(pvPoint))< 0.5 && sqrt(t->px()*t->px() + t->py()*t->py() ) > 0.4  && deltaR(t->eta(), t->phi(), LVTau.Eta(), LVTau.Phi()) < 1.){

            if (DEBUG) cout<<"Found track within Isolation"<<endl;

            std::vector<int>   iIsolationTrack_VertexWithSignalMuonIsValid;
            std::vector<float> iIsolationTrack_VertexWithSignalMuonChi2;
            std::vector<float> iIsolationTrack_VertexWithSignalMuonPosition;
            IsolationTracksIndices.push_back(iTrack);
            std::vector<float> iIsolation_Track_p4;
            iIsolation_Track_p4.push_back(sqrt(pow(t->p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
            iIsolation_Track_p4.push_back(t->px());
            iIsolation_Track_p4.push_back(t->py());
            iIsolation_Track_p4.push_back(t->pz());

            iIsolationTrack_charge.push_back(t->charge());
            iIsolationTrack_isHighPurity.push_back(t->quality(TrackBase::highPurity));

            IsolationTrack_p4.at(IsolationTrack_p4.size() - 1).push_back(iIsolation_Track_p4);

            iIsolationTrack_dxySV.push_back(t->dxy(TheSecondaryVertexPoint));
            iIsolationTrack_dzSV.push_back(t->dz(TheSecondaryVertexPoint));

            iIsolationTrack_dxyPV.push_back(t->dxy(pvPoint));
            iIsolationTrack_dzPV.push_back(t->dz(pvPoint));


            ClosestApproachInRPhi DocaMuon1, DocaMuon2, DocaMuon3;

            DocaMuon1.calculate(theB->build(t).initialFreeState(), iTransientTracks.at(sortedindices.at(0)).initialFreeState());
            DocaMuon2.calculate(theB->build(t).initialFreeState(), iTransientTracks.at(sortedindices.at(1)).initialFreeState());
            DocaMuon3.calculate(theB->build(t).initialFreeState(), iTransientTracks.at(sortedindices.at(2)).initialFreeState());


            if(DocaMuon1.status()){ iIsolationTrack_DocaMu1.push_back(DocaMuon1.distance());}
            else iIsolationTrack_DocaMu1.push_back(-1);

            if(DocaMuon2.status()){ iIsolationTrack_DocaMu2.push_back(DocaMuon2.distance());}
            else iIsolationTrack_DocaMu2.push_back(-1);

            if(DocaMuon3.status()){ iIsolationTrack_DocaMu3.push_back(DocaMuon3.distance());}
            else iIsolationTrack_DocaMu3.push_back(-1);
            iIsolationTracksCollection.push_back(TheB->build(t));

            vector<TransientTrack> IsolationTrack_Muon1;
            vector<TransientTrack> IsolationTrack_Muon2;
            vector<TransientTrack> IsolationTrack_Muon3;


            IsolationTrack_Muon1.push_back(iTransientTracks.at(0));IsolationTrack_Muon1.push_back(TheB->build(t));
            IsolationTrack_Muon2.push_back(iTransientTracks.at(1));IsolationTrack_Muon2.push_back(TheB->build(t));
            IsolationTrack_Muon3.push_back(iTransientTracks.at(2));IsolationTrack_Muon3.push_back(TheB->build(t));

            KalmanVertexFitter KVF_IsoTrack_Mu1(true),KVF_IsoTrack_Mu2(true),KVF_IsoTrack_Mu3(true); 
            TransientVertex V_IsoTrack_Mu1,V_IsoTrack_Mu2,V_IsoTrack_Mu3;
            bool FitTrackMu1Ok(true);
            bool FitTrackMu2Ok(true);
            bool FitTrackMu3Ok(true);

            try {
               V_IsoTrack_Mu1=KVF_IsoTrack_Mu1.vertex(IsolationTrack_Muon1);
            } catch (...) {
               FitTrackMu1Ok = false;
            }

            try {
               V_IsoTrack_Mu2=KVF_IsoTrack_Mu2.vertex(IsolationTrack_Muon2);
            } catch (...) {
               FitTrackMu2Ok = false;
            }

            try {
               V_IsoTrack_Mu3=KVF_IsoTrack_Mu3.vertex(IsolationTrack_Muon3);
            } catch (...) {
               FitTrackMu3Ok = false;
            }

            iIsolationTrack_VertexWithSignalMuonIsValid.push_back(V_IsoTrack_Mu1.isValid());
            iIsolationTrack_VertexWithSignalMuonIsValid.push_back(V_IsoTrack_Mu2.isValid());
            iIsolationTrack_VertexWithSignalMuonIsValid.push_back(V_IsoTrack_Mu3.isValid());
            if(FitTrackMu1Ok&&FitTrackMu2Ok&&FitTrackMu3Ok);// dummy; to make it compile
            if(V_IsoTrack_Mu1.isValid())
            {
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(V_IsoTrack_Mu1.totalChiSquared());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu1.position().x());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu1.position().y());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu1.position().z());
            }
            else 
            {
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(-1);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
            }

            if(V_IsoTrack_Mu2.isValid())
            {	    
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(V_IsoTrack_Mu2.totalChiSquared());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu2.position().x());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu2.position().y());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu2.position().z());
            }
            else 
            {
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(-1);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
            }

            if(V_IsoTrack_Mu3.isValid())
            {	    
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(V_IsoTrack_Mu3.totalChiSquared());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu3.position().x());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu3.position().y());
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(V_IsoTrack_Mu3.position().z());
            }
            else 
            {
               iIsolationTrack_VertexWithSignalMuonChi2.push_back(-1);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
               iIsolationTrack_VertexWithSignalMuonPosition.push_back(-99.);
            }
            IsolationTrack_VertexWithSignalMuonIsValid.at(IsolationTrack_VertexWithSignalMuonIsValid.size() - 1).push_back(iIsolationTrack_VertexWithSignalMuonIsValid);
            IsolationTrack_VertexWithSignalMuonChi2.at(IsolationTrack_VertexWithSignalMuonChi2.size() - 1).push_back(iIsolationTrack_VertexWithSignalMuonChi2);
            IsolationTrack_VertexWithSignalMuonPosition.at(IsolationTrack_VertexWithSignalMuonPosition.size() - 1).push_back(iIsolationTrack_VertexWithSignalMuonPosition);
         }



         //--------------------------- Isolation Branch

         double dz = abs(t->dz(TheSecondaryVertexPoint));
         double dxy = abs(t->dxy(TheSecondaryVertexPoint));
         double dca_fv = sqrt(dz*dz+dxy*dxy);
         double dr_tau = deltaR(t->eta(), t->phi(), LVTau.Eta(), LVTau.Phi());

         ////////////////////////////////////////////////
         // Below is the isolation defined by Jian
         // iso no. 1b - using pt_min, drtau_max of the 3 mu
         if(t->pt() > 0.33*minmuon_pt && dr_tau < 3.*maxmuon_dr && dca_fv<0.05 ) {
            sumptalltracks += t->pt();
            sumalltracks++; // iso 3b
            if(dca_fv<mindist)mindist=dca_fv; // iso 4b
         } 

         if(t->pt()<1.0) continue;  // was 1.2
         // iso no. 1
         if(dr_tau < 0.5 && dca_fv<0.05 ) {
            sumptalltracks05 += t->pt();
            sumalltracks05++; // iso 3
            if(dca_fv<mindca_iso05)mindca_iso05=dca_fv; // iso 4
         }

         if(dca_fv<mindca_iso05)mindist05=dca_fv; // iso 4
         if(dca_fv<0.05)sumalltracks_b++; // iso 3b


         TransientTrack trkiso = theB->build(t);
         ClosestApproachInRPhi cAppm1, cAppm2, cAppm3;
         cAppm1.calculate(trkiso.initialFreeState(), iTransientTracks.at(0).initialFreeState());
         cAppm2.calculate(trkiso.initialFreeState(), iTransientTracks.at(1).initialFreeState());
         cAppm3.calculate(trkiso.initialFreeState(), iTransientTracks.at(2).initialFreeState());
         if(!(cAppm1.status()&&cAppm2.status()&&cAppm3.status())) continue;


         // iso no. 2
         if(deltaR(t->eta(), t->phi(), LV1.Eta(), LV1.Phi()) < 0.3 && cAppm1.distance() < 0.1) {// && dz1 < .3) 
            N_trk_1++;
            pt_trk_1 += t->pt();
         }
         if(deltaR(t->eta(), t->phi(), LV2.Eta(), LV2.Phi()) < 0.3 && cAppm2.distance() < 0.1) {//&& dz2 < .3) 
            N_trk_2++;
            pt_trk_2 += t->pt();
         }
         if(deltaR(t->eta(), t->phi(), LV3.Eta(), LV3.Phi()) < 0.3 && cAppm3.distance() < 0.1) {//&& dz3 < .3) 
            N_trk_3++;
            pt_trk_3 += t->pt();
         }
         if( (deltaR(t->eta(), t->phi(), LV1.Eta(), LV1.Phi()) < 0.3 && cAppm1.distance() < 0.1 )
               ||(deltaR(t->eta(), t->phi(), LV2.Eta(), LV2.Phi()) < 0.3 && cAppm2.distance() < 0.1 )
               ||(deltaR(t->eta(), t->phi(), LV3.Eta(), LV3.Phi()) < 0.3 && cAppm3.distance() < 0.1 )
           ) N_trk_total++;



         double dz_primaryvertex=abs(t->dz(pvPoint));

         if(!(dz_primaryvertex < 1))continue;
         double dxy_primaryvertex = abs(t->dxy(pvPoint));
         if(dxy_primaryvertex>0.1) N_trk0p1++;
         if(dxy_primaryvertex>0.2) N_trk0p2++;
         if(dxy_primaryvertex>0.5) N_trk0p5++;
         if(dxy_primaryvertex>maxdxy) maxdxy = dxy_primaryvertex;
      }
      //********************************************************************************
      //Here reconstruct the vertices of all signal candidates with all isolation tracks 


      IsolationTrack_dxySV.push_back(iIsolationTrack_dxySV);
      IsolationTrack_dzSV.push_back(iIsolationTrack_dzSV);

      IsolationTrack_dxyPV.push_back(iIsolationTrack_dxyPV);
      IsolationTrack_dzPV.push_back(iIsolationTrack_dzPV);
      IsolationTrack_DocaMu1.push_back(iIsolationTrack_DocaMu1);
      IsolationTrack_DocaMu2.push_back(iIsolationTrack_DocaMu2);
      IsolationTrack_DocaMu3.push_back(iIsolationTrack_DocaMu3);
      IsolationTrack_charge.push_back(iIsolationTrack_charge);
      IsolationTrack_isHighPurity.push_back(iIsolationTrack_isHighPurity);

      relative_iso = sumptalltracks/LVTau.Pt();
      relative_iso05 = sumptalltracks05/LVTau.Pt();
      relative_mu1_iso = pt_trk_1/LV1.Pt(); relative_mu2_iso = pt_trk_2/LV2.Pt(); relative_mu3_iso = pt_trk_3/LV3.Pt();
      relative_maxiso = TMath::Max(relative_mu1_iso, TMath::Max(relative_mu2_iso,relative_mu3_iso ));

      std::vector<float> isolation1, isolation2, isolation3, isolation4;
      isolation1.push_back(relative_iso);
      isolation1.push_back(sumalltracks);
      isolation1.push_back(mindist);

      isolation2.push_back(relative_iso05);
      isolation2.push_back(sumalltracks05);
      isolation2.push_back(mindist);
      isolation2.push_back(mindist05);


      isolation3.push_back(relative_mu1_iso);
      isolation3.push_back(relative_mu2_iso);
      isolation3.push_back(relative_mu3_iso);
      isolation3.push_back(relative_maxiso);

      isolation4.push_back(N_trk_1);
      isolation4.push_back(N_trk_2);
      isolation4.push_back(N_trk_3);

      isolation4.push_back(N_trk0p1);
      isolation4.push_back(N_trk0p2);
      isolation4.push_back(N_trk0p5);
      isolation4.push_back(maxdxy);

      Vertex_Isolation1.push_back(isolation1);
      Vertex_Isolation2.push_back(isolation2);
      Vertex_Isolation3.push_back(isolation3);
      Vertex_Isolation4.push_back(isolation4);
      index++;
   }

   for(size_t isv = 0; isv < svertices->size(); isv++) {

      const reco::VertexCompositePtrCandidate& sv = (*svertices)[isv];
      SV_Track_P4.push_back(std::vector<std::vector<float> >());

      std::vector<float> iSV_pos;
      iSV_pos.push_back(sv.position().x());
      iSV_pos.push_back(sv.position().y());
      iSV_pos.push_back(sv.position().z());
      SV_pos.push_back(iSV_pos);
      SV_Mass.push_back(sv.p4().M());

      TMatrixTSym<float> sv_cov(LorentzVectorParticle::NVertex);

      math::Error<3>::type sv_Cov;
      sv.fillVertexCovariance(sv_Cov);

      for (int i = 0; i <LorentzVectorParticle::NVertex; i++){
         for (int j = 0; j < LorentzVectorParticle::NVertex; j++) {
            sv_cov(i, j) = sv_Cov(i, j);
            sv_cov(j, i) = sv_Cov(i, j);
         }
      }

      std::vector<float>  sv_covariance;
      for (int i = 0; i < LorentzVectorParticle::NVertex; i++) {
         for (int j = i; j < LorentzVectorParticle::NVertex; j++) {
            sv_covariance.push_back(sv_cov(i, j));
         }
      }
      SV_PosCovariance.push_back(sv_covariance);
      std::vector<int>    iSV_Trackcharge;

      size_t nDaughters = sv.numberOfDaughters();

      for (size_t itrk=0; itrk<nDaughters; ++itrk){ 
         std::vector<float>  iSV_TrackP4;

         const reco::Candidate & trk = *sv.daughter(itrk);

         iSV_TrackP4.push_back(sqrt(pow(trk.p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
         iSV_TrackP4.push_back(trk.px());
         iSV_TrackP4.push_back(trk.py());
         iSV_TrackP4.push_back(trk.pz());
         SV_Track_P4.at(SV_Track_P4.size()-1).push_back(iSV_TrackP4);
         iSV_Trackcharge.push_back(trk.charge());

      }
      SV_TrackCharge.push_back(iSV_Trackcharge);
   }

   return;
}
