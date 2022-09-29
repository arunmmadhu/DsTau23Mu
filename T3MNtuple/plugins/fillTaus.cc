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
        
           
        for(size_t i=0; i<cands->size(); ++i){
          if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;
          if(!(*cands)[i].bestTrack()) continue;
          
          unsigned int key = (*cands)[i].vertexRef().key();
          int quality = (*cands)[i].pvAssociationQuality();
        
        
          if(key!=0 || (quality!=pat::PackedCandidate::UsedInFitTight && quality!=pat::PackedCandidate::UsedInFitLoose)) continue;
          
          pvTracks.push_back(*((*cands)[i].bestTrack()));
        }
        
        
        
        
        // For the packedPFCandidate TRACKS from pfCandHandle, saved to pvertexTracks and allTracks
        //
        const edm::View<pat::PackedCandidate>* track_cands = tracksHandle.product();
        
        reco::TrackCollection pvertexTracks;
        reco::TrackCollection allTracks; // for taus (with possible SV) (testing now)
        
        
        for(size_t i=0; i<track_cands->size(); ++i)
          {
            if((*track_cands)[i].charge()==0 || (*track_cands)[i].vertexRef().isNull()) continue; // Make sure vertex collection isn't empty?
            if(!(*track_cands)[i].bestTrack()) continue; //if pointer doesn't exist, reject
            
            unsigned int key = (*track_cands)[i].vertexRef().key();
            int quality = (*track_cands)[i].pvAssociationQuality();
            
            // here I need to select "good" tracks
            // save them to all tracks
            // and if they belong to PV save them to pv tracks
            if (!(key!=0 ||
            (quality!=pat::PackedCandidate::UsedInFitTight
             && quality!=pat::PackedCandidate::UsedInFitLoose)))// key should be 0 and quality should be tight or loose?
                {
                  pvertexTracks.push_back(*((*track_cands)[i].bestTrack()));
                  // allTracks.push_back(*((*track_cands)[i].bestTrack())); // test for HelixLine Momentum is zero
                }
            
            // TODO: add requirement of "goodness"?
            allTracks.push_back(*((*track_cands)[i].bestTrack()));
            
            //cout << "Charged track with vertex filled. " << endl;
            
          }
          
          
          //Dealing with TauHandle
          for(unsigned iTau = 0; iTau < tauHandle->size(); iTau++){
            pat::TauRef tau(tauHandle,iTau);
            float valueAOD = tau->tauID("byIsolationMVArun2v1DBoldDMwLTraw");
            float valueMiniAOD = tau->tauID("byIsolationMVArun2v1DBoldDMwLTrawNew");//(*mvaIsoRaw)[tau];
            
            cout << "valueAOD: " << valueAOD << endl;
            cout << "valueMiniAOD: " << valueMiniAOD << endl;
            
            cout << "pT: " << tau->pt() << " eta: " << tau->eta() << " phi: " << tau->phi() << " energy: " << tau->energy() << endl;
            
            //if(valueAOD != valueMiniAOD) unmatchedTaus.push_back(tau);
          }

  }

