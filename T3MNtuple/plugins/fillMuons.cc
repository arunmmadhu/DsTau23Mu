#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
#include "DsTau23Mu/T3MNtuple/interface/TrackParticle.h"
#include "DsTau23Mu/T3MNtuple/interface/ParticleBuilder.h"

void T3MNtuple::fillMuons( const edm::Event& iEvent,
            const edm::EventSetup& iSetup,
            const Handle<std::vector<reco::Muon> >& muons,
            const Handle<TrackCollection>& trackCollection,
            const Handle<VertexCollection>& pvs)
      {

         unsigned int Muon_index = 0;
         unsigned int sel_muon_index = 0;
         ThreeMuons_index.resize(ThreeMuons_idx.size());
         TwoMuonsTrack_Muonsindex.resize(TwoMuonsTrack_idx.size());


         for (vector<reco::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon, Muon_index++) {
            reco::MuonRef RefMuon(muons, Muon_index);
            if(AcceptedMuon(RefMuon)){
               for(unsigned int iThreeMuons=0;  iThreeMuons < ThreeMuons_idx.size(); iThreeMuons++){
                  if(find(ThreeMuons_idx.at(iThreeMuons).begin(), ThreeMuons_idx.at(iThreeMuons).end(), Muon_index) !=  ThreeMuons_idx.at(iThreeMuons).end()){
                     ThreeMuons_index.at(iThreeMuons).push_back(sel_muon_index);
                  }
               }

               for(unsigned int iTwoMuons=0;  iTwoMuons < TwoMuonsTrack_idx.size(); iTwoMuons++){
                  if(TwoMuonsTrack_idx.at(iTwoMuons).at(0) == Muon_index || TwoMuonsTrack_idx.at(iTwoMuons).at(1) == Muon_index){
                     TwoMuonsTrack_Muonsindex.at(iTwoMuons).push_back(sel_muon_index);
                  }
               }

               std::vector<double> iMuon_Poca;
               iMuon_Poca.push_back(RefMuon->vx());
               iMuon_Poca.push_back(RefMuon->vy());
               iMuon_Poca.push_back(RefMuon->vz());
               Muon_Poca.push_back(iMuon_Poca);
               std::vector<double> iMuon_p4;
               iMuon_p4.push_back(RefMuon->p4().E());
               iMuon_p4.push_back(RefMuon->p4().Px());
               iMuon_p4.push_back(RefMuon->p4().Py());
               iMuon_p4.push_back(RefMuon->p4().Pz());
               Muon_p4.push_back(iMuon_p4);

               const MuonIsolation Iso03 = RefMuon->isolationR03();
               const MuonIsolation Iso05 = RefMuon->isolationR05();

               const MuonPFIsolation PFIso03 = RefMuon->pfIsolationR03();
               const MuonPFIsolation PFIso04 = RefMuon->pfIsolationR04();

               //Muon_matchedToKs.push_back(MuonMatchedKs(RefMuon, const std::vector<reco::VertexCompositePtrCandidate>& kshort);

               Muon_numberOfChambers.push_back(RefMuon->numberOfChambers());
               Muon_isGlobalMuon.push_back(RefMuon->isGlobalMuon());
               Muon_isPFMuon.push_back(RefMuon->isPFMuon());
               Muon_isRPCMuon.push_back(RefMuon->isRPCMuon());
               Muon_isStandAloneMuon.push_back(RefMuon->isStandAloneMuon());
               Muon_isTrackerMuon.push_back(RefMuon->isTrackerMuon());
               Muon_isCaloMuon.push_back(RefMuon->isCaloMuon());
               Muon_isQualityValid.push_back(RefMuon->isQualityValid());
               Muon_isTimeValid.push_back(RefMuon->isTimeValid());
               Muon_isIsolationValid.push_back(RefMuon->isIsolationValid());
               Muon_numberOfMatchedStations.push_back(RefMuon->numberOfMatchedStations());
               Muon_numberOfMatches.push_back(RefMuon->numberOfMatches(Muon::SegmentArbitration));
               Muon_charge.push_back(RefMuon->charge());

               Muon_expectedNnumberOfMatchedStations.push_back(RefMuon->expectedNnumberOfMatchedStations());


               const Vertex & VertexMuonID = (*pvs)[dump_pv_index_to_fill.at(0)];
               int idbit(0);
               if(muon::isLooseMuon(*RefMuon)) idbit |= 1 << 0;
               if(muon::isSoftMuon(*RefMuon,VertexMuonID)) idbit |= 1 << 1;
               if(muon::isMediumMuon(*RefMuon)) idbit |= 1 << 2;
               if(muon::isTightMuon(*RefMuon,VertexMuonID)) idbit |= 1 << 3;
               if(muon::isHighPtMuon(*RefMuon,VertexMuonID)) idbit |= 1 << 4;
               Muon_ID.push_back(idbit);
               //	       std::cout<<"muonsoftmva  "<< RefMuon->softMvaValue() << std::endl;  // only in Mini AOD

               int ssbit(0);
               if(RefMuon->passed(Muon::CutBasedIdLoose))ssbit|=1<<0;
               if(RefMuon->passed(Muon::CutBasedIdMedium))ssbit|=1<<1;
               if(RefMuon->passed(Muon::CutBasedIdMediumPrompt))ssbit|=1<<2;
               if(RefMuon->passed(Muon::CutBasedIdTight))ssbit|=1<<3;
               if(RefMuon->passed(Muon::CutBasedIdGlobalHighPt))ssbit|=1<<4;
               if(RefMuon->passed(Muon::CutBasedIdTrkHighPt))ssbit|=1<<5;
               if(RefMuon->passed(Muon::PFIsoVeryLoose))ssbit|=1<<6;
               if(RefMuon->passed(Muon::PFIsoLoose))ssbit|=1<<7;
               if(RefMuon->passed(Muon::PFIsoMedium))ssbit|=1<<8;
               if(RefMuon->passed(Muon::PFIsoTight))ssbit|=1<<9;
               if(RefMuon->passed(Muon::PFIsoVeryTight))ssbit|=1<<10;
               if(RefMuon->passed(Muon::TkIsoLoose))ssbit|=1<<11;
               if(RefMuon->passed(Muon::TkIsoTight))ssbit|=1<<12;
               if(RefMuon->passed(Muon::SoftCutBasedId))ssbit|=1<<13;
               if(RefMuon->passed(Muon::SoftMvaId))ssbit|=1<<14;
               if(RefMuon->passed(Muon::MvaLoose))ssbit|=1<<15;
               if(RefMuon->passed(Muon::MvaMedium))ssbit|=1<<16; 
               if(RefMuon->passed(Muon::MvaTight))ssbit|=1<<17;
               if(RefMuon->passed(Muon::MiniIsoLoose))ssbit|=1<<18;
               if(RefMuon->passed(Muon::MiniIsoMedium))ssbit|=1<<19;
               if(RefMuon->passed(Muon::MiniIsoTight))ssbit|=1<<20;
               if(RefMuon->passed(Muon::MiniIsoVeryTight))ssbit|=1<<21;


               Muon_StandardSelection.push_back(ssbit);
               /////////////////////////////////////////////////////////////
               //here following guide given in:
               //https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
               /////////////////////////////////////////////////////////////


               std::vector<float>  iMuon_TrackX;
               std::vector<float>  iMuon_TrackY;
               std::vector<float>  iMuon_dDxDz;
               std::vector<float>  iMuon_dDyDz;
               std::vector<float>  iMuon_dX;
               std::vector<float>  iMuon_dY;
               std::vector<float>  iMuon_pullX;
               std::vector<float>  iMuon_pullY;
               std::vector<float>  iMuon_pullDxDz;
               std::vector<float>  iMuon_pullDyDz;
               std::vector<float>  inumberOfSegments;


               for(int it = 1; it <= 8; ++it) { // loop over stations, 1,2,3,4 are DT, 5,6,7,8 are CSC
                  if(it <=4){ // we are in DT
                     iMuon_TrackX.push_back(RefMuon->trackX(it,1));
                     iMuon_TrackY.push_back(RefMuon->trackY(it,1));
                     iMuon_dDxDz.push_back(RefMuon->dDxDz(it,1));
                     iMuon_dDyDz.push_back(RefMuon->dDyDz(it,1));
                     iMuon_dX.push_back(RefMuon->dX(it,1));
                     iMuon_dY.push_back(RefMuon->dY(it,1));
                     iMuon_pullX.push_back(RefMuon->pullX(it,1));
                     iMuon_pullY.push_back(RefMuon->pullY(it,1));
                     iMuon_pullDxDz.push_back(RefMuon->pullDxDz(it,1));
                     iMuon_pullDyDz.push_back(RefMuon->pullDyDz(it,1));
                     inumberOfSegments.push_back(RefMuon->numberOfSegments(it,1));
                  }
                  if(it > 4){  // now in csc
                     iMuon_TrackX.push_back(RefMuon->trackX(it - 4,2));
                     iMuon_TrackY.push_back(RefMuon->trackY(it - 4,2));
                     iMuon_dDxDz.push_back(RefMuon->dDxDz(it - 4,2));
                     iMuon_dDyDz.push_back(RefMuon->dDyDz(it - 4,2));
                     iMuon_dX.push_back(RefMuon->dX(it - 4,2));
                     iMuon_dY.push_back(RefMuon->dY(it - 4,2));
                     iMuon_pullX.push_back(RefMuon->pullX(it - 4,2));
                     iMuon_pullY.push_back(RefMuon->pullY(it - 4,2));
                     iMuon_pullDxDz.push_back(RefMuon->pullDxDz(it - 4,2));
                     iMuon_pullDyDz.push_back(RefMuon->pullDyDz(it - 4,2));
                     inumberOfSegments.push_back(RefMuon->numberOfSegments(it - 4,2));
                  }
               }


               Muon_TrackX.push_back(iMuon_TrackX);
               Muon_TrackY.push_back(iMuon_TrackY);
               Muon_dDxDz.push_back(iMuon_dDxDz);
               Muon_dDyDz.push_back(iMuon_dDyDz);
               Muon_dX.push_back(iMuon_dX);
               Muon_dY.push_back(iMuon_dY);
               Muon_pullX.push_back(iMuon_pullX);
               Muon_pullY.push_back(iMuon_pullY);
               Muon_pullDxDz.push_back(iMuon_pullDxDz);
               Muon_pullDyDz.push_back(iMuon_pullDyDz);
               numberOfSegments.push_back(inumberOfSegments);

               std::vector<double> iMuon_outerTrack_p4;
               std::vector<double> iMuon_innerTrack_p4;
               if (RefMuon->isGlobalMuon()) {
                  Muon_hitPattern_numberOfValidMuonHits.push_back(RefMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
                  Muon_normChi2.push_back(RefMuon->globalTrack()->normalizedChi2());

                  iMuon_outerTrack_p4.push_back(RefMuon->outerTrack()->eta());
                  iMuon_outerTrack_p4.push_back(RefMuon->outerTrack()->phi());
                  Muon_prod_inner_outer_charge.push_back(RefMuon->outerTrack()->charge()*RefMuon->innerTrack()->charge());

                  Muon_outerTrack_normalizedChi2.push_back(RefMuon->outerTrack()->normalizedChi2());
                  Muon_outerTrack_muonStationsWithValidHits.push_back(RefMuon->outerTrack()->hitPattern().muonStationsWithValidHits());


                  unsigned int dt1(0),dt2(0),dt3(0),dt4(0);
                  unsigned int rpc1(0),rpc2(0),rpc3(0),rpc4(0);
                  unsigned int csc1(0),csc2(0),csc3(0),csc4(0);
                  double comb(0);
                  const reco::HitPattern &pattern = RefMuon->globalTrack()->hitPattern();
                  for (int i=0;i<pattern.numberOfAllHits(reco::HitPattern::TRACK_HITS);i++)
                  { 
                     uint32_t hit = pattern.getHitPattern(reco::HitPattern::TRACK_HITS,i);
                     if (pattern.validHitFilter(hit) != 1) {continue;}
                     if (pattern.getMuonStation(hit) == 1)
                     { 
                        if (pattern.muonDTHitFilter(hit))  dt1++;
                        if (pattern.muonRPCHitFilter(hit)) rpc1++;
                        if (pattern.muonCSCHitFilter(hit)) csc1++;
                     }
                     else if (pattern.getMuonStation(hit) == 2)
                     { 
                        if (pattern.muonDTHitFilter(hit))  dt2++;
                        if (pattern.muonRPCHitFilter(hit)) rpc2++;
                        if (pattern.muonCSCHitFilter(hit)) csc2++;
                     }
                     else if (pattern.getMuonStation(hit) == 3)
                     { 
                        if (pattern.muonDTHitFilter(hit))  dt3++;
                        if (pattern.muonRPCHitFilter(hit)) rpc3++;
                        if (pattern.muonCSCHitFilter(hit)) csc3++;
                     }
                     else if (pattern.getMuonStation(hit) == 4)
                     { 
                        if (pattern.muonDTHitFilter(hit))  dt4++;
                        if (pattern.muonRPCHitFilter(hit)) rpc4++;
                        if (pattern.muonCSCHitFilter(hit)) csc4++;
                     }    
                  }
                  comb = (dt1+dt2+dt3+dt4)/2. + (rpc1+rpc2+rpc3+rpc4);
                  csc1>6 ? comb+=6 : comb+=csc1;
                  csc2>6 ? comb+=6 : comb+=csc2;
                  csc3>6 ? comb+=6 : comb+=csc3;
                  csc4>6 ? comb+=6 : comb+=csc4;
                  Muon_vmuonhitcomb_reco.push_back(comb);
                  Muon_rpchits_reco.push_back(rpc1+rpc2+rpc3+rpc4);


               } else {
                  Muon_normChi2.push_back(0);
                  Muon_hitPattern_numberOfValidMuonHits.push_back(0);
                  Muon_prod_inner_outer_charge.push_back(0);
                  Muon_outerTrack_normalizedChi2.push_back(0);
                  Muon_outerTrack_muonStationsWithValidHits.push_back(0);
                  Muon_vmuonhitcomb_reco.push_back(0);
                  Muon_rpchits_reco.push_back(0);
               }

               if (RefMuon->isTrackerMuon() || RefMuon->isGlobalMuon()) {
                  Muon_numberofValidPixelHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits());
                  Muon_trackerLayersWithMeasurement.push_back(RefMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
                  Muon_innerTrack_validFraction.push_back(RefMuon->innerTrack()->validFraction());
                  Muon_innerTrack_pixelLayersWithMeasurement.push_back(RefMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement() );
                  Muon_innerTrack_numberOfValidTrackerHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfValidTrackerHits() );
                  Muon_innerTrack_numberOfLostTrackerHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS) );
                  Muon_innerTrack_numberOfLostTrackerInnerHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS) );
                  Muon_innerTrack_numberOfLostTrackerOuterHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS) );
                  Muon_innerTrack_normalizedChi2.push_back(RefMuon->innerTrack()->normalizedChi2() );

                  Muon_innerTrack_numberofValidHits.push_back(RefMuon->innerTrack()->numberOfValidHits());
                  Muon_hitPattern_pixelLayerwithMeas.push_back(RefMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());

                  Muon_innerTrack_quality.push_back(RefMuon->innerTrack()->quality(TrackBase::highPurity));
                  iMuon_innerTrack_p4.push_back(RefMuon->innerTrack()->p());
                  iMuon_innerTrack_p4.push_back(RefMuon->innerTrack()->eta());
                  iMuon_innerTrack_p4.push_back(RefMuon->innerTrack()->phi());
               } else {
                  Muon_trackerLayersWithMeasurement.push_back(0);
                  Muon_numberofValidPixelHits.push_back(0);
                  Muon_innerTrack_quality.push_back(0);
                  Muon_innerTrack_numberofValidHits.push_back(0);
                  Muon_hitPattern_pixelLayerwithMeas.push_back(0);
                  Muon_innerTrack_validFraction.push_back(0);
                  Muon_innerTrack_pixelLayersWithMeasurement.push_back(0);
                  Muon_innerTrack_numberOfValidTrackerHits.push_back(0);
                  Muon_innerTrack_numberOfLostTrackerHits.push_back(0);
                  Muon_innerTrack_numberOfLostTrackerInnerHits.push_back(0);
                  Muon_innerTrack_numberOfLostTrackerOuterHits.push_back(0);
                  Muon_innerTrack_normalizedChi2.push_back(0);
               }
               Muon_outerTrack_p4.push_back(iMuon_outerTrack_p4);
               Muon_innerTrack_p4.push_back(iMuon_innerTrack_p4);


               Muon_timeAtIpInOut.push_back(RefMuon->time().timeAtIpInOut);
               Muon_timeAtIpInOutErr.push_back(RefMuon->time().timeAtIpInOutErr);


               if (RefMuon->isIsolationValid()) {
                  Muon_emEt03.push_back(Iso03.emEt);
                  Muon_emVetoEt03.push_back(Iso03.emVetoEt);
                  Muon_hadEt03.push_back(Iso03.hadEt);
                  Muon_hadVetoEt03.push_back(Iso03.hadVetoEt);
                  Muon_nJets03.push_back(Iso03.nJets);
                  Muon_nTracks03.push_back(Iso03.nTracks);
                  Muon_sumPt03.push_back(Iso03.sumPt);
                  Muon_trackerVetoPt03.push_back(Iso03.trackerVetoPt);

                  Muon_emEt05.push_back(Iso05.emEt);
                  Muon_emVetoEt05.push_back(Iso05.emVetoEt);
                  Muon_hadEt05.push_back(Iso05.hadEt);
                  Muon_hadVetoEt05.push_back(Iso05.hadVetoEt);
                  Muon_nJets05.push_back(Iso05.nJets);
                  Muon_nTracks05.push_back(Iso05.nTracks);
                  Muon_sumPt05.push_back(Iso05.sumPt);
                  Muon_trackerVetoPt05.push_back(Iso05.trackerVetoPt);
               } else { // if isolation is not valid use -1 as default
                  Muon_emEt03.push_back(-1);
                  Muon_emVetoEt03.push_back(-1);
                  Muon_hadEt03.push_back(-1);
                  Muon_hadVetoEt03.push_back(-1);
                  Muon_nJets03.push_back(-1);
                  Muon_nTracks03.push_back(-1);
                  Muon_sumPt03.push_back(-1);
                  Muon_trackerVetoPt03.push_back(-1);

                  Muon_emEt05.push_back(-1);
                  Muon_emVetoEt05.push_back(-1);
                  Muon_hadEt05.push_back(-1);
                  Muon_hadVetoEt05.push_back(-1);
                  Muon_nJets05.push_back(-1);
                  Muon_nTracks05.push_back(-1);
                  Muon_sumPt05.push_back(-1);
                  Muon_trackerVetoPt05.push_back(-1);
               }


               //--- Fill1 PFMuonIsolation -----
               if (RefMuon->isPFIsolationValid()) {
                  Muon_sumChargedHadronPt03.push_back(PFIso03.sumChargedHadronPt);
                  Muon_sumChargedParticlePt03.push_back(PFIso03.sumChargedParticlePt);
                  Muon_sumNeutralHadronEt03.push_back(PFIso03.sumNeutralHadronEt);
                  Muon_sumNeutralHadronEtHighThreshold03.push_back(PFIso03.sumNeutralHadronEtHighThreshold);
                  Muon_sumPhotonEt03.push_back(PFIso03.sumPhotonEt);
                  Muon_sumPhotonEtHighThreshold03.push_back(PFIso03.sumPhotonEtHighThreshold);
                  Muon_sumPUPt03.push_back(PFIso03.sumPUPt);

                  Muon_sumChargedHadronPt04.push_back(PFIso04.sumChargedHadronPt);
                  Muon_sumChargedParticlePt04.push_back(PFIso04.sumChargedParticlePt);
                  Muon_sumNeutralHadronEt04.push_back(PFIso04.sumNeutralHadronEt);
                  Muon_sumNeutralHadronEtHighThreshold04.push_back(PFIso04.sumNeutralHadronEtHighThreshold);
                  Muon_sumPhotonEt04.push_back(PFIso04.sumPhotonEt);
                  Muon_sumPhotonEtHighThreshold04.push_back(PFIso04.sumPhotonEtHighThreshold);
                  Muon_sumPUPt04.push_back(PFIso04.sumPUPt);
               } else { // if isolation is not valid use -1 as default
                  Muon_sumChargedHadronPt03.push_back(-1);
                  Muon_sumChargedParticlePt03.push_back(-1);
                  Muon_sumNeutralHadronEt03.push_back(-1);
                  Muon_sumNeutralHadronEtHighThreshold03.push_back(-1);
                  Muon_sumPhotonEt03.push_back(-1);
                  Muon_sumPhotonEtHighThreshold03.push_back(-1);
                  Muon_sumPUPt03.push_back(-1);

                  Muon_sumChargedHadronPt04.push_back(-1);
                  Muon_sumChargedParticlePt04.push_back(-1);
                  Muon_sumNeutralHadronEt04.push_back(-1);
                  Muon_sumNeutralHadronEtHighThreshold04.push_back(-1);
                  Muon_sumPhotonEt04.push_back(-1);
                  Muon_sumPhotonEtHighThreshold04.push_back(-1);
                  Muon_sumPUPt04.push_back(-1);
               }


               ///////////////////////////////////// Muon Combined Quality /////////////////////////////////////////////////////////////////////////////////////
               //   find more about combined Muon quality in http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_4_4/doc/html/d4/d52/structreco_1_1MuonQuality.html
               Muon_combinedQuality_updatedSta.push_back(RefMuon->combinedQuality().updatedSta);
               Muon_combinedQuality_trkKink.push_back(RefMuon->combinedQuality().trkKink);
               Muon_combinedQuality_glbKink.push_back(RefMuon->combinedQuality().glbKink);
               Muon_combinedQuality_trkRelChi2.push_back(RefMuon->combinedQuality().trkRelChi2);
               Muon_combinedQuality_staRelChi2.push_back(RefMuon->combinedQuality().staRelChi2);
               Muon_combinedQuality_chi2LocalPosition.push_back(RefMuon->combinedQuality().chi2LocalPosition);
               Muon_combinedQuality_chi2LocalMomentum.push_back(RefMuon->combinedQuality().chi2LocalMomentum);
               Muon_combinedQuality_localDistance.push_back(RefMuon->combinedQuality().localDistance);
               Muon_combinedQuality_globalDeltaEtaPhi.push_back(RefMuon->combinedQuality().globalDeltaEtaPhi);
               Muon_combinedQuality_tightMatch.push_back(RefMuon->combinedQuality().tightMatch);
               Muon_combinedQuality_glbTrackProbability.push_back(RefMuon->combinedQuality().glbTrackProbability);

               Muon_calEnergy_em.push_back(RefMuon->calEnergy().em);
               Muon_calEnergy_emS9.push_back(RefMuon->calEnergy().emS9);
               Muon_calEnergy_emS25.push_back(RefMuon->calEnergy().emS25);
               Muon_calEnergy_had.push_back(RefMuon->calEnergy().had);
               Muon_calEnergy_hadS9.push_back(RefMuon->calEnergy().hadS9);

               Muon_segmentCompatibility.push_back(muon::segmentCompatibility(*RefMuon));
               Muon_caloCompatibility.push_back(muon::caloCompatibility(*RefMuon));

               Muon_ptErrOverPt.push_back(RefMuon->muonBestTrack()->ptError()/RefMuon->muonBestTrack()->pt());

               Muon_ptError.push_back(RefMuon->muonBestTrack()->ptError());
               Muon_phiError.push_back(RefMuon->muonBestTrack()->phiError());
               Muon_etaError.push_back(RefMuon->muonBestTrack()->etaError());

               Muon_isGoodMuon_TM2DCompatibility.push_back(muon::isGoodMuon(*RefMuon, muon::TM2DCompatibilityTight));
               Muon_isGoodMuon_TrackerMuonArbitrated.push_back(muon::isGoodMuon(*RefMuon,muon::TrackerMuonArbitrated));
               Muon_isGoodMuon_TMOneStationTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMOneStationTight));
               Muon_isGoodMuon_TMOneStationAngTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMOneStationAngTight));
               Muon_isGoodMuon_TMLastStationTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMLastStationTight));
               Muon_isGoodMuon_TMLastStationAngTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMLastStationAngTight));
               Muon_isGoodMuon_TMLastStationOptimizedLowPtTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMLastStationOptimizedLowPtTight));
               Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMLastStationOptimizedBarrelLowPtTight));


               reco::TrackRef Track = RefMuon->track();
               int ntp = Muon_par.size();
               Muon_par.push_back(std::vector<float>());
               Muon_cov.push_back(std::vector<float>());
               if (Track.isNonnull()) {
                  GlobalPoint pvpoint(Track->vx(), Track->vy(), Track->vz());
		  //                  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
                  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);
		  
		  const TransientTrackBuilder* transTrackBuilder = &iSetup.getData(ttkToken_);
                  reco::TransientTrack transTrk = transTrackBuilder->build(Track);

                  TrackParticle trackparticle = ParticleBuilder::CreateTrackParticle(transTrk, transTrackBuilder, pvpoint, true, true);
                  Muon_trackCharge.push_back(trackparticle.Charge());
                  Muon_pdgid.push_back(trackparticle.PDGID());
                  Muon_B.push_back(trackparticle.BField());
                  Muon_M.push_back(trackparticle.Mass());
		  //		  std::cout<<"Muon_M  "<< trackparticle.Mass() <<" B  "<<trackparticle.BField()<< " pdg  "<<trackparticle.PDGID() <<std::endl;
                  for (int i = 0; i < trackparticle.NParameters(); i++) {
                     Muon_par.at(ntp).push_back(trackparticle.Parameter(i));
		     //		     std::cout<<" par:  " << trackparticle.Parameter(i) <<std::endl;
		     //		     std::cout<<" muonpar size  "<< Muon_par.at(ntp).size() << std::endl;
                     for (int j = i; j < trackparticle.NParameters(); j++) {
		       Muon_cov.at(ntp).push_back(trackparticle.Covariance(i, j)); 
                     }
                  }
               } else {
                  Muon_trackCharge.push_back(-999);
                  Muon_pdgid.push_back(-999);
                  Muon_B.push_back(-999);
                  Muon_M.push_back(-999);
               }

               int match;
               getTrackMatch(trackCollection, Track, match);
               Muon_Track_idx.push_back(match);
               sel_muon_index++;
            }
         }

      }

void T3MNtuple::fillMuons( const edm::Event& iEvent,
            const edm::EventSetup& iSetup,
            const Handle<std::vector<pat::Muon> >& muons,
            const Handle<TrackCollection>& trackCollection,
            const Handle<VertexCollection>& pvs)
      {

         unsigned int Muon_index = 0;
         unsigned int sel_muon_index = 0;
         ThreeMuons_index.resize(ThreeMuons_idx.size());
         TwoMuonsTrack_Muonsindex.resize(TwoMuonsTrack_idx.size());


         for (vector<pat::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon, Muon_index++) {
            pat::MuonRef RefMuon(muons, Muon_index);
            if(AcceptedMuon(RefMuon)){
               for(unsigned int iThreeMuons=0;  iThreeMuons < ThreeMuons_idx.size(); iThreeMuons++){
                  if(find(ThreeMuons_idx.at(iThreeMuons).begin(), ThreeMuons_idx.at(iThreeMuons).end(), Muon_index) !=  ThreeMuons_idx.at(iThreeMuons).end()){
                     ThreeMuons_index.at(iThreeMuons).push_back(sel_muon_index);
                  }
               }

               for(unsigned int iTwoMuons=0;  iTwoMuons < TwoMuonsTrack_idx.size(); iTwoMuons++){
                  if(TwoMuonsTrack_idx.at(iTwoMuons).at(0) == Muon_index || TwoMuonsTrack_idx.at(iTwoMuons).at(1) == Muon_index){
                     TwoMuonsTrack_Muonsindex.at(iTwoMuons).push_back(sel_muon_index);
                  }
               }

               std::vector<double> iMuon_Poca;
               iMuon_Poca.push_back(RefMuon->vx());
               iMuon_Poca.push_back(RefMuon->vy());
               iMuon_Poca.push_back(RefMuon->vz());
               Muon_Poca.push_back(iMuon_Poca);
               std::vector<double> iMuon_p4;
               iMuon_p4.push_back(RefMuon->p4().E());
               iMuon_p4.push_back(RefMuon->p4().Px());
               iMuon_p4.push_back(RefMuon->p4().Py());
               iMuon_p4.push_back(RefMuon->p4().Pz());
               Muon_p4.push_back(iMuon_p4);

               const MuonIsolation Iso03 = RefMuon->isolationR03();
               const MuonIsolation Iso05 = RefMuon->isolationR05();

               const MuonPFIsolation PFIso03 = RefMuon->pfIsolationR03();
               const MuonPFIsolation PFIso04 = RefMuon->pfIsolationR04();

               //Muon_matchedToKs.push_back(MuonMatchedKs(RefMuon, const std::vector<reco::VertexCompositePtrCandidate>& kshort);

               Muon_numberOfChambers.push_back(RefMuon->numberOfChambers());
               Muon_isGlobalMuon.push_back(RefMuon->isGlobalMuon());
               Muon_isPFMuon.push_back(RefMuon->isPFMuon());
               Muon_isRPCMuon.push_back(RefMuon->isRPCMuon());
               Muon_isStandAloneMuon.push_back(RefMuon->isStandAloneMuon());
               Muon_isTrackerMuon.push_back(RefMuon->isTrackerMuon());
               Muon_isCaloMuon.push_back(RefMuon->isCaloMuon());
               Muon_isQualityValid.push_back(RefMuon->isQualityValid());
               Muon_isTimeValid.push_back(RefMuon->isTimeValid());
               Muon_isIsolationValid.push_back(RefMuon->isIsolationValid());
               Muon_numberOfMatchedStations.push_back(RefMuon->numberOfMatchedStations());
               Muon_numberOfMatches.push_back(RefMuon->numberOfMatches(Muon::SegmentArbitration));
               Muon_charge.push_back(RefMuon->charge());

               Muon_expectedNnumberOfMatchedStations.push_back(RefMuon->expectedNnumberOfMatchedStations());


               const Vertex & VertexMuonID = (*pvs)[dump_pv_index_to_fill.at(0)];
               int idbit(0);
               if(muon::isLooseMuon(*RefMuon)) idbit |= 1 << 0;
               if(muon::isSoftMuon(*RefMuon,VertexMuonID)) idbit |= 1 << 1;
               if(muon::isMediumMuon(*RefMuon)) idbit |= 1 << 2;
               if(muon::isTightMuon(*RefMuon,VertexMuonID)) idbit |= 1 << 3;
               if(muon::isHighPtMuon(*RefMuon,VertexMuonID)) idbit |= 1 << 4;
               Muon_ID.push_back(idbit);
               //	       std::cout<<"muonsoftmva  "<< RefMuon->softMvaValue() << std::endl;  // only in Mini AOD

               int ssbit(0);
               if(RefMuon->passed(Muon::CutBasedIdLoose))ssbit|=1<<0;
               if(RefMuon->passed(Muon::CutBasedIdMedium))ssbit|=1<<1;
               if(RefMuon->passed(Muon::CutBasedIdMediumPrompt))ssbit|=1<<2;
               if(RefMuon->passed(Muon::CutBasedIdTight))ssbit|=1<<3;
               if(RefMuon->passed(Muon::CutBasedIdGlobalHighPt))ssbit|=1<<4;
               if(RefMuon->passed(Muon::CutBasedIdTrkHighPt))ssbit|=1<<5;
               if(RefMuon->passed(Muon::PFIsoVeryLoose))ssbit|=1<<6;
               if(RefMuon->passed(Muon::PFIsoLoose))ssbit|=1<<7;
               if(RefMuon->passed(Muon::PFIsoMedium))ssbit|=1<<8;
               if(RefMuon->passed(Muon::PFIsoTight))ssbit|=1<<9;
               if(RefMuon->passed(Muon::PFIsoVeryTight))ssbit|=1<<10;
               if(RefMuon->passed(Muon::TkIsoLoose))ssbit|=1<<11;
               if(RefMuon->passed(Muon::TkIsoTight))ssbit|=1<<12;
               if(RefMuon->passed(Muon::SoftCutBasedId))ssbit|=1<<13;
               if(RefMuon->passed(Muon::SoftMvaId))ssbit|=1<<14;
               if(RefMuon->passed(Muon::MvaLoose))ssbit|=1<<15;
               if(RefMuon->passed(Muon::MvaMedium))ssbit|=1<<16; 
               if(RefMuon->passed(Muon::MvaTight))ssbit|=1<<17;
               if(RefMuon->passed(Muon::MiniIsoLoose))ssbit|=1<<18;
               if(RefMuon->passed(Muon::MiniIsoMedium))ssbit|=1<<19;
               if(RefMuon->passed(Muon::MiniIsoTight))ssbit|=1<<20;
               if(RefMuon->passed(Muon::MiniIsoVeryTight))ssbit|=1<<21;


               Muon_StandardSelection.push_back(ssbit);
               /////////////////////////////////////////////////////////////
               //here following guide given in:
               //https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
               /////////////////////////////////////////////////////////////


               std::vector<float>  iMuon_TrackX;
               std::vector<float>  iMuon_TrackY;
               std::vector<float>  iMuon_dDxDz;
               std::vector<float>  iMuon_dDyDz;
               std::vector<float>  iMuon_dX;
               std::vector<float>  iMuon_dY;
               std::vector<float>  iMuon_pullX;
               std::vector<float>  iMuon_pullY;
               std::vector<float>  iMuon_pullDxDz;
               std::vector<float>  iMuon_pullDyDz;
               std::vector<float>  inumberOfSegments;


               for(int it = 1; it <= 8; ++it) { // loop over stations, 1,2,3,4 are DT, 5,6,7,8 are CSC
                  if(it <=4){ // we are in DT
                     iMuon_TrackX.push_back(RefMuon->trackX(it,1));
                     iMuon_TrackY.push_back(RefMuon->trackY(it,1));
                     iMuon_dDxDz.push_back(RefMuon->dDxDz(it,1));
                     iMuon_dDyDz.push_back(RefMuon->dDyDz(it,1));
                     iMuon_dX.push_back(RefMuon->dX(it,1));
                     iMuon_dY.push_back(RefMuon->dY(it,1));
                     iMuon_pullX.push_back(RefMuon->pullX(it,1));
                     iMuon_pullY.push_back(RefMuon->pullY(it,1));
                     iMuon_pullDxDz.push_back(RefMuon->pullDxDz(it,1));
                     iMuon_pullDyDz.push_back(RefMuon->pullDyDz(it,1));
                     inumberOfSegments.push_back(RefMuon->numberOfSegments(it,1));
                  }
                  if(it > 4){  // now in csc
                     iMuon_TrackX.push_back(RefMuon->trackX(it - 4,2));
                     iMuon_TrackY.push_back(RefMuon->trackY(it - 4,2));
                     iMuon_dDxDz.push_back(RefMuon->dDxDz(it - 4,2));
                     iMuon_dDyDz.push_back(RefMuon->dDyDz(it - 4,2));
                     iMuon_dX.push_back(RefMuon->dX(it - 4,2));
                     iMuon_dY.push_back(RefMuon->dY(it - 4,2));
                     iMuon_pullX.push_back(RefMuon->pullX(it - 4,2));
                     iMuon_pullY.push_back(RefMuon->pullY(it - 4,2));
                     iMuon_pullDxDz.push_back(RefMuon->pullDxDz(it - 4,2));
                     iMuon_pullDyDz.push_back(RefMuon->pullDyDz(it - 4,2));
                     inumberOfSegments.push_back(RefMuon->numberOfSegments(it - 4,2));
                  }
               }


               Muon_TrackX.push_back(iMuon_TrackX);
               Muon_TrackY.push_back(iMuon_TrackY);
               Muon_dDxDz.push_back(iMuon_dDxDz);
               Muon_dDyDz.push_back(iMuon_dDyDz);
               Muon_dX.push_back(iMuon_dX);
               Muon_dY.push_back(iMuon_dY);
               Muon_pullX.push_back(iMuon_pullX);
               Muon_pullY.push_back(iMuon_pullY);
               Muon_pullDxDz.push_back(iMuon_pullDxDz);
               Muon_pullDyDz.push_back(iMuon_pullDyDz);
               numberOfSegments.push_back(inumberOfSegments);

               std::vector<double> iMuon_outerTrack_p4;
               std::vector<double> iMuon_innerTrack_p4;
               if (RefMuon->isGlobalMuon()) {
                  Muon_hitPattern_numberOfValidMuonHits.push_back(RefMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
                  Muon_normChi2.push_back(RefMuon->globalTrack()->normalizedChi2());

                  iMuon_outerTrack_p4.push_back(RefMuon->outerTrack()->eta());
                  iMuon_outerTrack_p4.push_back(RefMuon->outerTrack()->phi());
                  Muon_prod_inner_outer_charge.push_back(RefMuon->outerTrack()->charge()*RefMuon->innerTrack()->charge());

                  Muon_outerTrack_normalizedChi2.push_back(RefMuon->outerTrack()->normalizedChi2());
                  Muon_outerTrack_muonStationsWithValidHits.push_back(RefMuon->outerTrack()->hitPattern().muonStationsWithValidHits());


                  unsigned int dt1(0),dt2(0),dt3(0),dt4(0);
                  unsigned int rpc1(0),rpc2(0),rpc3(0),rpc4(0);
                  unsigned int csc1(0),csc2(0),csc3(0),csc4(0);
                  double comb(0);
                  const reco::HitPattern &pattern = RefMuon->globalTrack()->hitPattern();
                  for (int i=0;i<pattern.numberOfAllHits(reco::HitPattern::TRACK_HITS);i++)
                  { 
                     uint32_t hit = pattern.getHitPattern(reco::HitPattern::TRACK_HITS,i);
                     if (pattern.validHitFilter(hit) != 1) {continue;}
                     if (pattern.getMuonStation(hit) == 1)
                     { 
                        if (pattern.muonDTHitFilter(hit))  dt1++;
                        if (pattern.muonRPCHitFilter(hit)) rpc1++;
                        if (pattern.muonCSCHitFilter(hit)) csc1++;
                     }
                     else if (pattern.getMuonStation(hit) == 2)
                     { 
                        if (pattern.muonDTHitFilter(hit))  dt2++;
                        if (pattern.muonRPCHitFilter(hit)) rpc2++;
                        if (pattern.muonCSCHitFilter(hit)) csc2++;
                     }
                     else if (pattern.getMuonStation(hit) == 3)
                     { 
                        if (pattern.muonDTHitFilter(hit))  dt3++;
                        if (pattern.muonRPCHitFilter(hit)) rpc3++;
                        if (pattern.muonCSCHitFilter(hit)) csc3++;
                     }
                     else if (pattern.getMuonStation(hit) == 4)
                     { 
                        if (pattern.muonDTHitFilter(hit))  dt4++;
                        if (pattern.muonRPCHitFilter(hit)) rpc4++;
                        if (pattern.muonCSCHitFilter(hit)) csc4++;
                     }    
                  }
                  comb = (dt1+dt2+dt3+dt4)/2. + (rpc1+rpc2+rpc3+rpc4);
                  csc1>6 ? comb+=6 : comb+=csc1;
                  csc2>6 ? comb+=6 : comb+=csc2;
                  csc3>6 ? comb+=6 : comb+=csc3;
                  csc4>6 ? comb+=6 : comb+=csc4;
                  Muon_vmuonhitcomb_reco.push_back(comb);
                  Muon_rpchits_reco.push_back(rpc1+rpc2+rpc3+rpc4);


               } else {
                  Muon_normChi2.push_back(0);
                  Muon_hitPattern_numberOfValidMuonHits.push_back(0);
                  Muon_prod_inner_outer_charge.push_back(0);
                  Muon_outerTrack_normalizedChi2.push_back(0);
                  Muon_outerTrack_muonStationsWithValidHits.push_back(0);
                  Muon_vmuonhitcomb_reco.push_back(0);
                  Muon_rpchits_reco.push_back(0);
               }

               if (RefMuon->isTrackerMuon() || RefMuon->isGlobalMuon()) {
                  Muon_numberofValidPixelHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits());
                  Muon_trackerLayersWithMeasurement.push_back(RefMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
                  Muon_innerTrack_validFraction.push_back(RefMuon->innerTrack()->validFraction());
                  Muon_innerTrack_pixelLayersWithMeasurement.push_back(RefMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement() );
                  Muon_innerTrack_numberOfValidTrackerHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfValidTrackerHits() );
                  Muon_innerTrack_numberOfLostTrackerHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS) );
                  Muon_innerTrack_numberOfLostTrackerInnerHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS) );
                  Muon_innerTrack_numberOfLostTrackerOuterHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS) );
                  Muon_innerTrack_normalizedChi2.push_back(RefMuon->innerTrack()->normalizedChi2() );

                  Muon_innerTrack_numberofValidHits.push_back(RefMuon->innerTrack()->numberOfValidHits());
                  Muon_hitPattern_pixelLayerwithMeas.push_back(RefMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());

                  Muon_innerTrack_quality.push_back(RefMuon->innerTrack()->quality(TrackBase::highPurity));
                  iMuon_innerTrack_p4.push_back(RefMuon->innerTrack()->p());
                  iMuon_innerTrack_p4.push_back(RefMuon->innerTrack()->eta());
                  iMuon_innerTrack_p4.push_back(RefMuon->innerTrack()->phi());
               } else {
                  Muon_trackerLayersWithMeasurement.push_back(0);
                  Muon_numberofValidPixelHits.push_back(0);
                  Muon_innerTrack_quality.push_back(0);
                  Muon_innerTrack_numberofValidHits.push_back(0);
                  Muon_hitPattern_pixelLayerwithMeas.push_back(0);
                  Muon_innerTrack_validFraction.push_back(0);
                  Muon_innerTrack_pixelLayersWithMeasurement.push_back(0);
                  Muon_innerTrack_numberOfValidTrackerHits.push_back(0);
                  Muon_innerTrack_numberOfLostTrackerHits.push_back(0);
                  Muon_innerTrack_numberOfLostTrackerInnerHits.push_back(0);
                  Muon_innerTrack_numberOfLostTrackerOuterHits.push_back(0);
                  Muon_innerTrack_normalizedChi2.push_back(0);
               }
               Muon_outerTrack_p4.push_back(iMuon_outerTrack_p4);
               Muon_innerTrack_p4.push_back(iMuon_innerTrack_p4);


               Muon_timeAtIpInOut.push_back(RefMuon->time().timeAtIpInOut);
               Muon_timeAtIpInOutErr.push_back(RefMuon->time().timeAtIpInOutErr);


               if (RefMuon->isIsolationValid()) {
                  Muon_emEt03.push_back(Iso03.emEt);
                  Muon_emVetoEt03.push_back(Iso03.emVetoEt);
                  Muon_hadEt03.push_back(Iso03.hadEt);
                  Muon_hadVetoEt03.push_back(Iso03.hadVetoEt);
                  Muon_nJets03.push_back(Iso03.nJets);
                  Muon_nTracks03.push_back(Iso03.nTracks);
                  Muon_sumPt03.push_back(Iso03.sumPt);
                  Muon_trackerVetoPt03.push_back(Iso03.trackerVetoPt);

                  Muon_emEt05.push_back(Iso05.emEt);
                  Muon_emVetoEt05.push_back(Iso05.emVetoEt);
                  Muon_hadEt05.push_back(Iso05.hadEt);
                  Muon_hadVetoEt05.push_back(Iso05.hadVetoEt);
                  Muon_nJets05.push_back(Iso05.nJets);
                  Muon_nTracks05.push_back(Iso05.nTracks);
                  Muon_sumPt05.push_back(Iso05.sumPt);
                  Muon_trackerVetoPt05.push_back(Iso05.trackerVetoPt);
               } else { // if isolation is not valid use -1 as default
                  Muon_emEt03.push_back(-1);
                  Muon_emVetoEt03.push_back(-1);
                  Muon_hadEt03.push_back(-1);
                  Muon_hadVetoEt03.push_back(-1);
                  Muon_nJets03.push_back(-1);
                  Muon_nTracks03.push_back(-1);
                  Muon_sumPt03.push_back(-1);
                  Muon_trackerVetoPt03.push_back(-1);

                  Muon_emEt05.push_back(-1);
                  Muon_emVetoEt05.push_back(-1);
                  Muon_hadEt05.push_back(-1);
                  Muon_hadVetoEt05.push_back(-1);
                  Muon_nJets05.push_back(-1);
                  Muon_nTracks05.push_back(-1);
                  Muon_sumPt05.push_back(-1);
                  Muon_trackerVetoPt05.push_back(-1);
               }


               //--- Fill1 PFMuonIsolation -----
               if (RefMuon->isPFIsolationValid()) {
                  Muon_sumChargedHadronPt03.push_back(PFIso03.sumChargedHadronPt);
                  Muon_sumChargedParticlePt03.push_back(PFIso03.sumChargedParticlePt);
                  Muon_sumNeutralHadronEt03.push_back(PFIso03.sumNeutralHadronEt);
                  Muon_sumNeutralHadronEtHighThreshold03.push_back(PFIso03.sumNeutralHadronEtHighThreshold);
                  Muon_sumPhotonEt03.push_back(PFIso03.sumPhotonEt);
                  Muon_sumPhotonEtHighThreshold03.push_back(PFIso03.sumPhotonEtHighThreshold);
                  Muon_sumPUPt03.push_back(PFIso03.sumPUPt);

                  Muon_sumChargedHadronPt04.push_back(PFIso04.sumChargedHadronPt);
                  Muon_sumChargedParticlePt04.push_back(PFIso04.sumChargedParticlePt);
                  Muon_sumNeutralHadronEt04.push_back(PFIso04.sumNeutralHadronEt);
                  Muon_sumNeutralHadronEtHighThreshold04.push_back(PFIso04.sumNeutralHadronEtHighThreshold);
                  Muon_sumPhotonEt04.push_back(PFIso04.sumPhotonEt);
                  Muon_sumPhotonEtHighThreshold04.push_back(PFIso04.sumPhotonEtHighThreshold);
                  Muon_sumPUPt04.push_back(PFIso04.sumPUPt);
               } else { // if isolation is not valid use -1 as default
                  Muon_sumChargedHadronPt03.push_back(-1);
                  Muon_sumChargedParticlePt03.push_back(-1);
                  Muon_sumNeutralHadronEt03.push_back(-1);
                  Muon_sumNeutralHadronEtHighThreshold03.push_back(-1);
                  Muon_sumPhotonEt03.push_back(-1);
                  Muon_sumPhotonEtHighThreshold03.push_back(-1);
                  Muon_sumPUPt03.push_back(-1);

                  Muon_sumChargedHadronPt04.push_back(-1);
                  Muon_sumChargedParticlePt04.push_back(-1);
                  Muon_sumNeutralHadronEt04.push_back(-1);
                  Muon_sumNeutralHadronEtHighThreshold04.push_back(-1);
                  Muon_sumPhotonEt04.push_back(-1);
                  Muon_sumPhotonEtHighThreshold04.push_back(-1);
                  Muon_sumPUPt04.push_back(-1);
               }


               ///////////////////////////////////// Muon Combined Quality /////////////////////////////////////////////////////////////////////////////////////
               //   find more about combined Muon quality in http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_4_4/doc/html/d4/d52/structreco_1_1MuonQuality.html
               Muon_combinedQuality_updatedSta.push_back(RefMuon->combinedQuality().updatedSta);
               Muon_combinedQuality_trkKink.push_back(RefMuon->combinedQuality().trkKink);
               Muon_combinedQuality_glbKink.push_back(RefMuon->combinedQuality().glbKink);
               Muon_combinedQuality_trkRelChi2.push_back(RefMuon->combinedQuality().trkRelChi2);
               Muon_combinedQuality_staRelChi2.push_back(RefMuon->combinedQuality().staRelChi2);
               Muon_combinedQuality_chi2LocalPosition.push_back(RefMuon->combinedQuality().chi2LocalPosition);
               Muon_combinedQuality_chi2LocalMomentum.push_back(RefMuon->combinedQuality().chi2LocalMomentum);
               Muon_combinedQuality_localDistance.push_back(RefMuon->combinedQuality().localDistance);
               Muon_combinedQuality_globalDeltaEtaPhi.push_back(RefMuon->combinedQuality().globalDeltaEtaPhi);
               Muon_combinedQuality_tightMatch.push_back(RefMuon->combinedQuality().tightMatch);
               Muon_combinedQuality_glbTrackProbability.push_back(RefMuon->combinedQuality().glbTrackProbability);

               Muon_calEnergy_em.push_back(RefMuon->calEnergy().em);
               Muon_calEnergy_emS9.push_back(RefMuon->calEnergy().emS9);
               Muon_calEnergy_emS25.push_back(RefMuon->calEnergy().emS25);
               Muon_calEnergy_had.push_back(RefMuon->calEnergy().had);
               Muon_calEnergy_hadS9.push_back(RefMuon->calEnergy().hadS9);

               Muon_segmentCompatibility.push_back(muon::segmentCompatibility(*RefMuon));
               Muon_caloCompatibility.push_back(muon::caloCompatibility(*RefMuon));

               Muon_ptErrOverPt.push_back(RefMuon->muonBestTrack()->ptError()/RefMuon->muonBestTrack()->pt());

               Muon_ptError.push_back(RefMuon->muonBestTrack()->ptError());
               Muon_phiError.push_back(RefMuon->muonBestTrack()->phiError());
               Muon_etaError.push_back(RefMuon->muonBestTrack()->etaError());

               Muon_isGoodMuon_TM2DCompatibility.push_back(muon::isGoodMuon(*RefMuon, muon::TM2DCompatibilityTight));
               Muon_isGoodMuon_TrackerMuonArbitrated.push_back(muon::isGoodMuon(*RefMuon,muon::TrackerMuonArbitrated));
               Muon_isGoodMuon_TMOneStationTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMOneStationTight));
               Muon_isGoodMuon_TMOneStationAngTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMOneStationAngTight));
               Muon_isGoodMuon_TMLastStationTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMLastStationTight));
               Muon_isGoodMuon_TMLastStationAngTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMLastStationAngTight));
               Muon_isGoodMuon_TMLastStationOptimizedLowPtTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMLastStationOptimizedLowPtTight));
               Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight.push_back(muon::isGoodMuon(*RefMuon,muon::TMLastStationOptimizedBarrelLowPtTight));


               reco::TrackRef Track = RefMuon->track();
               int ntp = Muon_par.size();
               Muon_par.push_back(std::vector<float>());
               Muon_cov.push_back(std::vector<float>());
               if (Track.isNonnull()) {
                  GlobalPoint pvpoint(Track->vx(), Track->vy(), Track->vz());
		  //                  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
		  //                  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);
		  const TransientTrackBuilder* transTrackBuilder = &iSetup.getData(ttkToken_);
                  reco::TransientTrack transTrk = transTrackBuilder->build(Track);
                  TrackParticle trackparticle = ParticleBuilder::CreateTrackParticle(transTrk, transTrackBuilder, pvpoint, true, true);
                  Muon_trackCharge.push_back(trackparticle.Charge());
                  Muon_pdgid.push_back(trackparticle.PDGID());
                  Muon_B.push_back(trackparticle.BField());
                  Muon_M.push_back(trackparticle.Mass());
                  for (int i = 0; i < trackparticle.NParameters(); i++) {
                     Muon_par.at(ntp).push_back(trackparticle.Parameter(i));
                     for (int j = i; j < trackparticle.NParameters(); j++) {
		       Muon_cov.at(ntp).push_back(trackparticle.Covariance(i, j)); // comment out to keep sizee low, this is unused variable
                     }
                  }
               } else {
                  Muon_trackCharge.push_back(-999);
                  Muon_pdgid.push_back(-999);
                  Muon_B.push_back(-999);
                  Muon_M.push_back(-999);
               }

               int match;
               getTrackMatch(trackCollection, Track, match);
               Muon_Track_idx.push_back(match);
               sel_muon_index++;
            }
         }

      }
