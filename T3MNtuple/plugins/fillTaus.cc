#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
#include "DsTau23Mu/T3MNtuple/interface/TrackParticle.h"
#include "DsTau23Mu/T3MNtuple/interface/ParticleBuilder.h"

void T3MNtuple::fillTaus(const edm::Event& iEvent,
            const edm::EventSetup& iSetup,
            const Handle<TrackCollection>& trackCollection,
            const Handle<VertexCollection>& pvs,
            const Handle<BeamSpot>& beamSpotHandle,
            const Handle<pat::TauCollection>& tauHandle,
            const Handle<vector<Vertex> >&  vertexs,
            const Handle<edm::View<pat::PackedCandidate> >& pfCandHandle,
            const Handle<edm::View<pat::PackedCandidate> >& tracksHandle)
      {

        // PF Tau Discriminators: tauIntDiscrims_ and tauFloatDiscrims_
        //
        tauIntDiscrims_ = 
        {
          "decayModeFinding", // it is decayModeFindingOldDMs
          "decayModeFindingNewDMs",
          
          "byLooseCombinedIsolationDeltaBetaCorr3Hits",
          "byMediumCombinedIsolationDeltaBetaCorr3Hits",
          "byTightCombinedIsolationDeltaBetaCorr3Hits",
          
          "byVLooseIsolationMVArun2v1DBoldDMwLT",
          "byLooseIsolationMVArun2v1DBoldDMwLT",
          "byMediumIsolationMVArun2v1DBoldDMwLT",
          "byTightIsolationMVArun2v1DBoldDMwLT",
          "byVTightIsolationMVArun2v1DBoldDMwLT",
        
        
          "byVLooseIsolationMVArun2v1DBnewDMwLT",    
          "byLooseIsolationMVArun2v1DBnewDMwLT",
          "byMediumIsolationMVArun2v1DBnewDMwLT",
          "byTightIsolationMVArun2v1DBnewDMwLT",
          "byVTightIsolationMVArun2v1DBnewDMwLT",
        
        
          "byLooseIsolationMVArun2v1DBdR03oldDMwLT",
          "byMediumIsolationMVArun2v1DBdR03oldDMwLT",
          "byTightIsolationMVArun2v1DBdR03oldDMwLT",
          "byVTightIsolationMVArun2v1DBdR03oldDMwLT",
          
          "byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",
          "byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",
          "byTightCombinedIsolationDeltaBetaCorr3HitsdR03",
        
        
          "againstElectronMVA5category",
          
          "byLooseIsolationMVA3newDMwLT",
          "byLooseIsolationMVA3oldDMwLT",
          "byLoosePileupWeightedIsolation3Hits",
          "byMediumIsolationMVA3newDMwLT",
          "byMediumIsolationMVA3oldDMwLT",
          "byMediumPileupWeightedIsolation3Hits",
          "byTightIsolationMVA3newDMwLT",
          "byTightIsolationMVA3oldDMwLT",
          "byTightPileupWeightedIsolation3Hits",
          
          "byVLooseIsolationMVA3newDMwLT",
          "byVTightIsolationMVA3newDMwLT",
          "byVVTightIsolationMVA3newDMwLT",
        
        
          "byVLooseIsolationMVA3oldDMwLT",
          "byVTightIsolationMVA3oldDMwLT",
          "byVVTightIsolationMVA3oldDMwLT",
        
        
          "againstMuonLoose3",
          "againstMuonTight3",
        
        
          "againstElectronVLooseMVA6",
          "againstElectronLooseMVA6",
          "againstElectronMediumMVA6",
          "againstElectronTightMVA6",
          "againstElectronVTightMVA6"
        };
        
        
        tauFloatDiscrims_ =
        {
          "byCombinedIsolationDeltaBetaCorrRaw3Hits",
          "byIsolationMVArun2v1DBoldDMwLTraw",
          "byIsolationMVA3oldDMwoLTraw",
          "byIsolationMVA3oldDMwLTraw",
          "byIsolationMVA3newDMwoLTraw",
          "againstElectronMVA5raw",
          "byPhotonPtSumOutsideSignalCone",
          "byPileupWeightedIsolationRaw3Hits",
          "footprintCorrection",
          "neutralIsoPtSumWeight",
          "photonPtSumOutsideSignalCone",
          "byIsolationMVA3newDMwLTraw",
          "chargedIsoPtSum",
          "neutralIsoPtSum",
          "puCorrPtSum",
        };
        
        
        
        
        
        // For the packedPFCandidate TRACKS from pfCandHandle, saved to pvTracks
        //
        const edm::View<pat::PackedCandidate>* cands = pfCandHandle.product();
        TLorentzVector aTrack;
        reco::TrackCollection pvTracks;
        reco::TrackCollection allTracks;
        
           
        for(size_t i=0; i<cands->size(); ++i){
          if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;
          if(!(*cands)[i].bestTrack()) continue;
          
          unsigned int key = (*cands)[i].vertexRef().key();
          int quality = (*cands)[i].pvAssociationQuality();
        
	  allTracks.push_back(*((*cands)[i].bestTrack()));
	  // key == 0 means the tracks that are assigned to the First PV in the collection, i.e. key = vertex index
          if(key!=0 || (quality!=pat::PackedCandidate::UsedInFitTight && quality!=pat::PackedCandidate::UsedInFitLoose)) continue;
	  pvTracks.push_back(*((*cands)[i].bestTrack()));
        }
        
        
        
	const Vertex & Highest_pT_vertex = (*pvs)[0]; //  simply take the first in the collection;
	std::vector<float> iHighestPt_PrimaryVertex_Pos;
	if(Highest_pT_vertex.isValid())
	  {
	    iHighestPt_PrimaryVertex_Pos.push_back(Highest_pT_vertex.position().x());
	    iHighestPt_PrimaryVertex_Pos.push_back(Highest_pT_vertex.position().y());
	    iHighestPt_PrimaryVertex_Pos.push_back(Highest_pT_vertex.position().z());
	  }
	Vertex_HighestPt_PrimaryVertex.push_back(iHighestPt_PrimaryVertex_Pos);

	TMatrixTSym<double> hpTcov(LorentzVectorParticle::NVertex);
	math::Error<3>::type hpTCov;
	Highest_pT_vertex.fill(hpTCov);

	for (int i = 0; i <LorentzVectorParticle::NVertex; i++)
	  {
	  for (int j = 0; j < LorentzVectorParticle::NVertex; j++) 
	    {
	      hpTcov(i, j) = hpTCov(i, j);
	      hpTcov(j, i) = hpTCov(i, j);
	    }
	  }
	std::vector<double>  hpT_cov;     
	for (int i = 0; i < LorentzVectorParticle::NVertex; i++) 
	  {
	    for (int j = i; j < LorentzVectorParticle::NVertex; j++) 
	    {
	      hpT_cov.push_back(hpTcov(i, j));
	    }
	  }
	
	Vertex_HighestPt_PrimaryVertex_covariance.push_back(hpT_cov);
	

          //Dealing with TauHandle
          for(unsigned iTau = 0; iTau < tauHandle->size(); iTau++){
            pat::TauRef tau(tauHandle,iTau);
            float valueAOD = tau->tauID("byIsolationMVArun2v1DBoldDMwLTraw");
	    //            float valueMiniAOD = tau->tauID("byIsolationMVArun2v1DBoldDMwLTrawNew");//(*mvaIsoRaw)[tau];
            float valueMiniAOD = tau->tauID("byLooseDeepTau2017v2p1VSmu");//(*mvaIsoRaw)[tau];
            
	    // ---------------  fill momentum and ID variables --------------
	    std::vector<float> iTau_p4;
	    iTau_p4.push_back(tau->energy());
	    iTau_p4.push_back(tau->px());
	    iTau_p4.push_back(tau->py());
	    iTau_p4.push_back(tau->pz());
	    Tau_p4.push_back(iTau_p4);

	    Tau_charge.push_back(tau->charge());
	    Tau_DecayMode.push_back(tau->decayMode());
	    Tau_DecayModeFinding.push_back(tau->tauID("decayModeFinding"));

	    Tau_byLooseDeepTau2017v2p1VSe.push_back(tau->tauID("byLooseDeepTau2017v2p1VSe"));
	    Tau_byMediumDeepTau2017v2p1VSe.push_back(tau->tauID("byMediumDeepTau2017v2p1VSe"));
	    Tau_byTightDeepTau2017v2p1VSe.push_back(tau->tauID("byTightDeepTau2017v2p1VSe"));

	    Tau_byLooseDeepTau2017v2p1VSmu.push_back(tau->tauID("byLooseDeepTau2017v2p1VSmu"));
	    Tau_byMediumDeepTau2017v2p1VSmu.push_back(tau->tauID("byMediumDeepTau2017v2p1VSmu"));
	    Tau_byTightDeepTau2017v2p1VSmu.push_back(tau->tauID("byTightDeepTau2017v2p1VSmu"));

	    Tau_byLooseDeepTau2017v2p1VSjet.push_back(tau->tauID("byLooseDeepTau2017v2p1VSjet"));
	    Tau_byMediumDeepTau2017v2p1VSjet.push_back(tau->tauID("byMediumDeepTau2017v2p1VSjet"));
	    Tau_byTightDeepTau2017v2p1VSjet.push_back(tau->tauID("byTightDeepTau2017v2p1VSjet"));

	    std::cout<<"  charged iso cands size "<< tau->signalChargedHadrCands().size() <<"   neutral    " << tau->signalNeutrHadrCands().size() <<"  photons  " << tau ->signalGammaCands().size() <<"   et le decay mode:   " <<tau->decayMode() <<"    " <<tau->leadChargedHadrCand()->p4().eta() <<std::endl;


	    std::vector<double > PFTauTrackLV;    
	    edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
	    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);

	    if(tau->decayMode() == 0  or tau->decayMode() == 1 )  //    either tau -> pi nu or tau->pi pi0 nu
	      {


		float PFTauTrack_deltaR=999.;
		reco::CandidatePtrVector  chargedCandidate = tau->signalChargedHadrCands();

		double deltaR(999.); 
		reco::Track RefToTauTrack;

		//allTracks
		//		for(auto iter: pvTracks)
		for(auto iter: allTracks)
		  {
		    float iDr = sqrt(pow(iter.eta() - tau->leadChargedHadrCand()->p4().eta(),2) +
				     pow(iter.phi() - tau->leadChargedHadrCand()->p4().phi(),2));
		    if( iDr< deltaR)
		      {
			deltaR = iDr;
			RefToTauTrack = iter;
		      }
		  }
	      
	    
	    
		Tau_PFTauTrackLV.push_back(tau->leadChargedHadrCand()->p4().e());    
		Tau_PFTauTrackLV.push_back(tau->leadChargedHadrCand()->p4().px());    
		Tau_PFTauTrackLV.push_back(tau->leadChargedHadrCand()->p4().py());    
		Tau_PFTauTrackLV.push_back(tau->leadChargedHadrCand()->p4().pz());    

		const reco::Track *TauTrack  =  &RefToTauTrack;
		
		GlobalPoint pvpoint(TauTrack->vx(), TauTrack->vy(), TauTrack->vz());
		reco::TransientTrack transTrk = transTrackBuilder->build(TauTrack);
		TrackParticle tautrackparticle = ParticleBuilder::CreateTrackParticle(transTrk, transTrackBuilder, pvpoint, true, true);

	


		//	std::cout<<"  is the right track found ??  "<< PFTauTrack_deltaR << std::endl;
		if(deltaR< 0.01)  //if track is matched; just arbitrary value
		  {
		    Tau_Track_Charge=tautrackparticle.Charge();
		    Tau_Track_pdgid=tautrackparticle.PDGID();
		    Tau_Track_B=tautrackparticle.BField();
		    Tau_Track_M=tautrackparticle.Mass();
		    
		    for (int i = 0; i < tautrackparticle.NParameters(); i++) 
		      {
			Tau_Track_par.push_back(tautrackparticle.Parameter(i));
			for (int j = i; j <tautrackparticle.NParameters(); j++) 
			  {
			    Tau_Track_cov.push_back(tautrackparticle.Covariance(i, j));
			  }
		      }
		  }else{
		  Tau_Track_Charge=-999;
		  Tau_Track_pdgid=-999;
		  Tau_Track_B=-999;
		  Tau_Track_M=-999;
		}
	      }
	  
	    if (tau->decayMode() == 10 or tau->decayMode() == 11) 
	      {
		
		std::vector<reco::TransientTrack> transTrk;
		TransientVertex transVtx;
		// const reco::PFCandidateRefVector cands = l.signalChargedHadrCands();
		reco::CandidatePtrVector signalCandidates = tau->signalChargedHadrCands();//signalCands();
     
		//		edm::Handle<reco::BeamSpot> beamSpotHandle;
		//		iEvent.getByToken(beamSpotTag, beamSpot);
 

		std::vector<reco::TransientTrack> transTracks;  
		//   find  tracks belonging to tau decay
		for (reco::CandidatePtrVector::const_iterator itr = signalCandidates.begin(); itr != signalCandidates.end(); ++itr) 
		  {
		    double deltaR(999.); 
		    reco::Track closestTrack;
		    for(auto iter: allTracks) 
		      {
			if( sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2))  < deltaR)
			  {
			    deltaR = sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2));
			    closestTrack = iter;
			  }
		      }
		    if(closestTrack.pt()!=0)transTracks.push_back(transTrackBuilder->build(closestTrack));
		  }
		std::cout<<" trans tracks:   "<< transTracks.size() << std::endl;

		bool fitOk = false;  
		if(transTracks.size() >= 2 ) {
		  KalmanVertexFitter kvf;
		  try {
		    transVtx = kvf.vertex(transTracks);
		    fitOk = true; 
		  } catch (...) {
		    fitOk = false; 
		    std::cout<<"Vtx fit failed!"<<std::endl;
		  }
		}


		fitOk = fitOk && transVtx.isValid() && fabs(transVtx.position().x())<1 && fabs(transVtx.position().y())<1;
		std::cout<<" is fit OK ??  "<<   fitOk << << <<std::endl;
    
		if(fitOk) {


		  Tau_SVPos.push_back(transVtx.position().x());
		  Tau_SVPos.push_back(transVtx.position().y());
		  Tau_SVPos.push_back(transVtx.position().z());



		  reco::Vertex secondaryVertex = transVtx;
 
		  //		  Tau_SVChi2.push_back();
		  //		  Tau_SV_Status_Chi2.push_back(secondaryVertex.chi2());

		  TMatrixTSym<double> svcov(3);
		  math::Error<3>::type svCov;
		  secondaryVertex.fill(svCov);
		  for (int i = 0; i <3; i++){
		    for (int j = 0; j < 3; j++) {
		      svcov(i, j) = svCov(i, j);
		      svcov(j, i) = svCov(i, j);
		      // cout<<"  svcov  "<<svcov(i,j)<<endl;
		    }
		  }
		  for (int i = 0; i < 3; i++) {
		    for (int j = i; j < 3; j++) {
		      Tau_SVCov.push_back(svcov(i, j));
		    }
		  }

		}  // if(fitOk)


	      } //  if (tau->decayMode() == 10 or tau->decayMode() == 11) 




	  }
      }



