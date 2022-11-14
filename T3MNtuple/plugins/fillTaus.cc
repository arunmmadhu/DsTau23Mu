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
          
	  /*          "byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",
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
          "byVVTightIsolationMVA3oldDMwLT", */
                
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
          "byIsolationMVArun2v1DBoldDMwLTraw",
	  //          "byIsolationMVA3oldDMwoLTraw",
	  //          "byIsolationMVA3oldDMwLTraw",
	  //          "byIsolationMVA3newDMwoLTraw",
	  //          "againstElectronMVA5raw",
          "byPhotonPtSumOutsideSignalCone",
	  //          "byPileupWeightedIsolationRaw3Hits",
	  //          "footprintCorrection",
          "neutralIsoPtSumWeight",
          "photonPtSumOutsideSignalCone",
	  //          "byIsolationMVA3newDMwLTraw",
          "chargedIsoPtSum",
          "neutralIsoPtSum",
          "puCorrPtSum",
        };
        
        
        // For the packedPFCandidate TRACKS from pfCandHandle, saved to pvTracks

        const edm::View<pat::PackedCandidate>* cands = pfCandHandle.product();
        TLorentzVector aTrack;
        reco::TrackCollection pvTracks;
        reco::TrackCollection allTracks;
        
           
        for(size_t i=0; i<cands->size(); ++i){
          if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;
          if(!(*cands)[i].bestTrack()) continue;
          
          unsigned int key = (*cands)[i].vertexRef().key();
          int quality = (*cands)[i].pvAssociationQuality();
        
	  allTracks.push_back(*((*cands)[i].bestTrack())); 	  // key == 0 means the tracks that are assigned to the First PV in the collection, i.e. key = vertex index
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
	
	
	for(unsigned iTau = 0; iTau < tauHandle->size(); iTau++){
	  pat::TauRef tau(tauHandle,iTau);
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
	  Tau_NewDecayModeFinding.push_back(tau->tauID("decayModeFindingNewDMs"));
	  //--------------------------   Selected Tau ID's 
	  Tau_byLooseDeepTau2017v2p1VSe.push_back(tau->tauID("byLooseDeepTau2017v2p1VSe"));
	  Tau_byMediumDeepTau2017v2p1VSe.push_back(tau->tauID("byMediumDeepTau2017v2p1VSe"));
	  Tau_byTightDeepTau2017v2p1VSe.push_back(tau->tauID("byTightDeepTau2017v2p1VSe"));
	  
	  Tau_byLooseDeepTau2017v2p1VSmu.push_back(tau->tauID("byLooseDeepTau2017v2p1VSmu"));
	  Tau_byMediumDeepTau2017v2p1VSmu.push_back(tau->tauID("byMediumDeepTau2017v2p1VSmu"));
	  Tau_byTightDeepTau2017v2p1VSmu.push_back(tau->tauID("byTightDeepTau2017v2p1VSmu"));

	  Tau_byLooseDeepTau2017v2p1VSjet.push_back(tau->tauID("byLooseDeepTau2017v2p1VSjet"));
	  Tau_byMediumDeepTau2017v2p1VSjet.push_back(tau->tauID("byMediumDeepTau2017v2p1VSjet"));
	  Tau_byTightDeepTau2017v2p1VSjet.push_back(tau->tauID("byTightDeepTau2017v2p1VSjet"));


	  Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
	  Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
	  Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));

	  Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));


	  //  fill the rest of float ID's 
	  std::vector<float> iTau_FloatDiscriminants;
	  for (unsigned int i = 0; i < tauFloatDiscrims_.size(); i++)
	    {
	      string ID = tauFloatDiscrims_.at(i);
	      float discr(-999.);
	      if(tau->isTauIDAvailable(ID.c_str())) discr = tau->tauID(ID.c_str());
	      iTau_FloatDiscriminants.push_back(discr);
	    }

	  // fill all userints
	  std::vector<int> iTau_IntDiscriminants;
	  for (unsigned int i = 0; i < tauIntDiscrims_.size(); i++)
	    {
	      string ID = tauIntDiscrims_.at(i);
	      float discr(-999.);
	      if(tau->isTauIDAvailable(ID.c_str())) discr = tau->tauID(ID.c_str());
	      iTau_IntDiscriminants.push_back(discr);

	    }


	  Tau_FloatDiscriminants.push_back(iTau_FloatDiscriminants);	    
	  Tau_IntDiscriminants.push_back(iTau_IntDiscriminants);


	  std::vector<double > PFTauTrackLV;    
	  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
	  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);


	  int NTau = Tau_a1_lvp.size();
	  Tau_SVPos.push_back(std::vector<float>());
	  Tau_SVCov.push_back(std::vector<float>());
	  Tau_a1_lvp.push_back(std::vector<float>());
	  Tau_a1_cov.push_back(std::vector<float>());
	  Tau_a1_charge.push_back(std::vector<int>());
	  Tau_a1_pdgid.push_back(std::vector<int>());
	  Tau_a1_B.push_back(std::vector<float>());
	  Tau_a1_M.push_back(std::vector<float>());
	  


	  Tau_PFTauTrack_p4.push_back(std::vector<float>());
	  Tau_Track_par.push_back(std::vector<float>());
	  Tau_Track_cov.push_back(std::vector<float>());
	  Tau_Track_Charge.push_back(std::vector<int>());
	  Tau_Track_pdgid.push_back(std::vector<int>());
	  Tau_Track_B.push_back(std::vector<float>());
	  Tau_Track_M.push_back(std::vector<float>());


	  if(tau->decayMode() == 0  or tau->decayMode() == 1 )  //    either tau -> pi nu or tau->pi pi0 nu
	    {
	      reco::CandidatePtrVector  chargedCandidate = tau->signalChargedHadrCands();
	      double deltaR(999.); 
	      reco::Track RefToTauTrack;
	      //  for(auto iter: pvTracks)
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
	      
	      Tau_PFTauTrack_p4.at(NTau).push_back(tau->leadChargedHadrCand()->p4().e());    
	      Tau_PFTauTrack_p4.at(NTau).push_back(tau->leadChargedHadrCand()->p4().px());    
	      Tau_PFTauTrack_p4.at(NTau).push_back(tau->leadChargedHadrCand()->p4().py());    
	      Tau_PFTauTrack_p4.at(NTau).push_back(tau->leadChargedHadrCand()->p4().pz());    
	      
	      const reco::Track *TauTrack  =  &RefToTauTrack;
	      
	      GlobalPoint pvpoint(TauTrack->vx(), TauTrack->vy(), TauTrack->vz());
	      reco::TransientTrack transTrk = transTrackBuilder->build(TauTrack);
	      TrackParticle tautrackparticle = ParticleBuilder::CreateTrackParticle(transTrk, transTrackBuilder, pvpoint, true, true);
	      
	      
	      if(deltaR< 0.01)  //if track is matched; just arbitrary value
		{
		  Tau_Track_Charge.at(NTau).push_back(tautrackparticle.Charge());
		  Tau_Track_pdgid.at(NTau).push_back(tautrackparticle.PDGID());
		  Tau_Track_B.at(NTau).push_back(tautrackparticle.BField());
		  Tau_Track_M.at(NTau).push_back(tautrackparticle.Mass());
		  
		  for (int i = 0; i < tautrackparticle.NParameters(); i++) 
		    {
		      Tau_Track_par.at(NTau).push_back(tautrackparticle.Parameter(i));
		      for (int j = i; j <tautrackparticle.NParameters(); j++) 
			{
			  Tau_Track_cov.at(NTau).push_back(tautrackparticle.Covariance(i, j));
			}
		    }
		}
	    }
	  
	  
	    if (tau->decayMode() == 10 or tau->decayMode() == 11) 
	      {
		
		std::vector<reco::TransientTrack> transTrk;
		TransientVertex transVtx;
		reco::CandidatePtrVector signalCandidates = tau->signalChargedHadrCands();//signalCands();
		std::vector<reco::TransientTrack> transTracks;  
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

		if(fitOk) {
		  Tau_SVPos.at(NTau).push_back(transVtx.position().x());
		  Tau_SVPos.at(NTau).push_back(transVtx.position().y());
		  Tau_SVPos.at(NTau).push_back(transVtx.position().z());

		  reco::Vertex secondaryVertex = transVtx;
		  TMatrixTSym<double> svcov(3);
		  math::Error<3>::type svCov;
		  secondaryVertex.fill(svCov);
		  for (int i = 0; i <3; i++)
		    {
		    for (int j = 0; j < 3; j++) 
		      {
			svcov(i, j) = svCov(i, j);
			svcov(j, i) = svCov(i, j);
		      }
		    }
		  for (int i = 0; i < 3; i++) 
		    {
		    for (int j = i; j < 3; j++) 
		      {
			Tau_SVCov.at(NTau).push_back(svcov(i, j));
		      }
		    }
		  

		  LorentzVectorParticle a1;
		  GlobalPoint sv(secondaryVertex.position().x(), secondaryVertex.position().y(), secondaryVertex.position().z());
		  KinematicParticleFactoryFromTransientTrack kinFactory;
		  float piMassSigma(sqrt(pow(10., -12.))), piChi(0.0), piNdf(0.0);
		  std::vector<RefCountedKinematicParticle> pions;
		  for (unsigned int i = 0; i <transTracks.size(); i++)
		    pions.push_back(kinFactory.particle(transTracks.at(i), PDGInfo::pi_mass(), piChi, piNdf, sv, piMassSigma));
      
		  KinematicParticleVertexFitter kpvFitter;
		  RefCountedKinematicTree jpTree = kpvFitter.fit(pions);
		  if(jpTree->isValid()){
		    jpTree->movePointerToTheTop();
		    const KinematicParameters parameters = jpTree->currentParticle()->currentState().kinematicParameters();
		    AlgebraicSymMatrix77 cov = jpTree->currentParticle()->currentState().kinematicParametersError().matrix();
		    TMatrixT<double> a1_par(LorentzVectorParticle::NLorentzandVertexPar, 1);
		    TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
		    for (int i = 0; i < LorentzVectorParticle::NLorentzandVertexPar; i++) 
		      {
			a1_par(i, 0) = parameters(i);
			for (int j = 0; j < LorentzVectorParticle::NLorentzandVertexPar; j++) 
			  {
			    a1_cov(i, j) = cov(i, j);
			  }
		      }
		    a1 = LorentzVectorParticle(a1_par, a1_cov, abs(PDGInfo::a_1_plus) * tau->charge(), tau->charge(), transTrackBuilder->field()->inInverseGeV(sv).z());
		    Tau_a1_charge.at(NTau).push_back(a1.Charge());
		    Tau_a1_pdgid.at(NTau).push_back(a1.PDGID());
		    Tau_a1_B.at(NTau).push_back(a1.BField());
		    Tau_a1_M.at(NTau).push_back(a1.Mass());
		    for (int i = 0; i < a1.NParameters(); i++) 
		      {
			Tau_a1_lvp.at(NTau).push_back(a1.Parameter(i));
			for (int j = i; j < a1.NParameters(); j++) 
			  {
			    Tau_a1_cov.at(NTau).push_back(a1.Covariance(i, j));
			  }
		      }
		  }
		  
		}  // if(fitOk)
		
	      } //  if (tau->decayMode() == 10 or tau->decayMode() == 11) 
	    
	  }
      }



