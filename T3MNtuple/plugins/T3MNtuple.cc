#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
//
// constants, enums and typedefs
//

//
// static data member definitions
//
double T3MNtuple::MuonPtCut_(-1.);
double T3MNtuple::MuonEtaCut_(999);
//
// constructors and destructor
//
T3MNtuple::T3MNtuple(const edm::ParameterSet& iConfig):
  muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  //badmuonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("badmuons"))),
  //jetToken_(consumes<reco::JetCollection>(iConfig.getParameter<edm::InputTag>("ak4PFJetsCHS"))),
  btagCvsBToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("btagsCvsB"))),
  btagCSVToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("btagsCSV"))),
  btagMVAToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("btagsMVA"))),
  vtxToken_(consumes<VertexCollection>(iConfig.getParameter<InputTag>("pvs"))),
  svToken_(consumes<VertexCollection>(iConfig.getParameter<InputTag>("svs"))),
  trackToken_(consumes<TrackCollection>(iConfig.getParameter<InputTag>("trks"))),
  triggerToken_(consumes<TriggerResults>(iConfig.getParameter<InputTag>("triggerBitsH"))),
  trigeventToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<InputTag>("triggerSummary"))),
  algToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<InputTag>("AlgInputTag"))),
  //level1Token_(consumes<L1GlobalTriggerReadoutRecord>(iConfig.getParameter<InputTag>("gtRecards"))),
  bsToken_(consumes<BeamSpot>(iConfig.getParameter<InputTag>("beamSpotHandle"))),
  puToken_(consumes<vector<PileupSummaryInfo> >(iConfig.getParameter<InputTag>("pileupSummary"))),
  genToken_(consumes<GenParticleCollection>(iConfig.getParameter<InputTag>("genParticles")))
  //BadGlbMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadGlbMuonFilter")))
  //refitter_(iConfig)
{

  gtUtil_ = new L1TGlobalUtil(iConfig, consumesCollector(), *this, algInputTag_, algInputTag_);
  //ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("PFCaloCompatibility");
  //muonCaloCompatibility_.configure( parameters );

  //ParameterSet parameters2 = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  //ConsumesCollector iC = consumesCollector();
  //parameters_.loadParameters( parameters2, iC );

  //trackAssociator_.useDefaultPropagator();

  doMC_ = iConfig.getParameter<bool>("doMC");
  wideSB_ = iConfig.getParameter<bool>("wideSB");
  do2mu_ = iConfig.getParameter<bool>("do2mu");
  passhlt_ = iConfig.getParameter<bool>("passhlt");
  mid_ = iConfig.getParameter<int>("mid");
  doTracks_ = iConfig.getParameter<bool>("doTracks");
  doMuons_ = iConfig.getParameter<bool>("doMuons");
  do3mutuple_ = iConfig.getParameter<bool>("do3mutuple");
  doL1_ = iConfig.getParameter<bool>("doL1");

  MuonPtCut_ = iConfig.getParameter<double>("MuonPtCut"); //default: 3.0
  MuonEtaCut_ = iConfig.getParameter<double>("MuonEtaCut"); //default: 2.5


}


T3MNtuple::~T3MNtuple()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//


bool T3MNtuple::isGoodTrack(const Track &track) {
  if(!(track.pt()>1.)){
    if(!(abs(track.eta())<2.4)){
      if(!(track.hitPattern().trackerLayersWithMeasurement()>5)){
	if(!(track.hitPattern().pixelLayersWithMeasurement()>1)) return false;
      }
    }
  }
  return true;
}



// ------------ method called for each event  ------------
void
T3MNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
  //Handle<bool> ifilterbadGlbMuon;
  //iEvent.getByToken(BadGlbMuonFilterToken_, ifilterbadGlbMuon);
  //filterbadGlbMuon = *ifilterbadGlbMuon;
  //cout<<filterbadGlbMuon<<endl;
  cnt_++;
  ClearEvent();


  fillEventInfo(iEvent, iSetup);
  if(doTracks_)
    fillTracks(iEvent, iSetup);
  if(doMuons_)
    fillMuons(iEvent, iSetup);
  if (doMC_)
    fillMCTruth(iEvent, iSetup);
  if (doL1_)
    fillL1(iEvent, iSetup);



  h_step->Fill(0);

  BeamSpot bs;
  Handle<BeamSpot> beamSpotHandle;
  iEvent.getByToken(bsToken_, beamSpotHandle);
  bs = *beamSpotHandle;

  //cout<<bs.x0()<<"  "<<bs.y0()<<"  "<<bs.z0()<<endl;
  //cout<<beamSpotHandle->position().x()<<"  "<<beamSpotHandle->position().y()<<"  "<<beamSpotHandle->position().z()<<endl;

  Handle<VertexCollection> pvs;
  iEvent.getByToken(vtxToken_ , pvs);
  n_vtx = pvs->size();

  Handle<VertexCollection> svs;
  iEvent.getByToken(svToken_ , svs);

  //cout<<(*pvs)[0].position().X()<<" "<<(*pvs)[0].position().Y()<<"  "<<(*pvs)[0].position().Z()<<endl;

  Handle<MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  if(muons->size() < (do2mu_?2:3))return;
  //Handle<ValueMap<MuonShower> > muonShowerInformationValueMapH_;
  //iEvent.getByToken(MuonShowerInformationValueMapToken_, muonShowerInformationValueMapH_);

  //Handle<MuonCollection> badmuons;
  //iEvent.getByToken(badmuonToken_, badmuons);
  //cout<<badmuons->size()<<endl;


  h_step->Fill(1);

  Handle<TrackCollection> trks;
  iEvent.getByToken(trackToken_, trks);

  //Handle<JetCollection> jets;
  //iEvents.getByToken(jetToken_, jets);
  Handle<JetTagCollection> btagsCvsB;
  iEvent.getByToken(btagCvsBToken_, btagsCvsB);
  Handle<JetTagCollection> btagsCSV;
  iEvent.getByToken(btagCSVToken_, btagsCSV);
  Handle<JetTagCollection> btagsMVA;
  iEvent.getByToken(btagMVAToken_, btagsMVA);

  Handle<TriggerResults> triggerBitsH;
  iEvent.getByToken(triggerToken_, triggerBitsH);

  Handle<trigger::TriggerEvent> triggerSummary;
  iEvent.getByToken(trigeventToken_, triggerSummary);

  Handle<BXVector<GlobalAlgBlk>> alg;
  iEvent.getByToken(algToken_,alg);





  //////////////////////////////
  // HLT

  hlt_doublemu4_lmnrt = 0;
  hlt_doublemu3_tau3mu = 0;

  const TriggerNames &triggerNames = iEvent.triggerNames( *triggerBitsH );
  for (size_t i_hlt = 0; i_hlt != triggerBitsH->size(); ++i_hlt)
    {
      string hltName = triggerNames.triggerName(i_hlt);
      //cout<<hltName<<"\t"<<triggerBitsH->accept(i_hlt)<<endl;
      if(!(hltName.find("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v") == string::npos)){
        if( triggerBitsH->wasrun(i_hlt) && !triggerBitsH->error(i_hlt) && triggerBitsH->accept(i_hlt )) hlt_doublemu4_lmnrt = 1;
      }
      if(!(hltName.find("HLT_DoubleMu3_Trk_Tau3mu_") == string::npos)){
        if( triggerBitsH->wasrun(i_hlt) && !triggerBitsH->error(i_hlt) && triggerBitsH->accept(i_hlt )) hlt_doublemu3_tau3mu = 1;
      }
    }
  if(passhlt_ && hlt_doublemu3_tau3mu<0.1) return;
  h_step->Fill(2);


    //////////////////////
    // b, c origin ?
  gen_flavor=1; nmu_mom=0;
  double gen_pv = 0;
  size_t ndsgen = 0;

  if(doMC_) {

    Handle<GenParticleCollection> genParticles;
    iEvent.getByToken(genToken_, genParticles);
    const GenParticle & gen1st = (*genParticles)[2];
    gen_pv = gen1st.vz();

    for(size_t i = 2; i < genParticles->size(); ++ i) {

      const GenParticle & p = (*genParticles)[i];
      //const Candidate * mom = p.mother();
      //if(i==2)cout<<setprecision(3)<< setiosflags(ios::showpoint)
      //  <<"id: "<<p.pdgId()<<"\tstatus: "<<p.status()
      //  <<"\tpt: "<<p.pt()<<"\teta: "<<p.eta()
      //  <<"\tphi: "<<p.phi()<<"\tmass: "<<p.mass()
      //  <<"\tvx: "<<p.vx()<<"\tvy "<<p.vy()
      //  <<"\tvz: "<<p.vz()<<"\tMOTHER id: "<<mom->pdgId()
      //  <<"\tMOTHER pt:"<<mom->pt()<<endl<<endl;

      if(abs(p.pdgId())==4 && gen_flavor==1){gen_flavor=4;}
      if(abs(p.pdgId())==5 && gen_flavor==1){gen_flavor=5;}

      if(abs(p.pdgId())==431) ndsgen++; // Ds
      if(abs(p.pdgId())!=13) continue; // mu
      if(abs(p.mother()->pdgId())!=mid_) continue; // phi (norm. channel), or tau (signal channel)
      nmu_mom++;
    }
    //if(ndsgen>1)return;
    //if(nmu_mom>3)return; // why ???
  }

  //h_step->Fill(2);

  ///////////////////////////
  // start to loop over muons
  n_reco = 0;
  pdgid_reco[0] = 0; pdgid_reco[1] = 0; pdgid_reco[2] = 0;
  momid_reco[0] = 0; momid_reco[1] = 0; momid_reco[2] = 0;
  vxy_reco[0] = 0; vxy_reco[1] = 0; vxy_reco[2] = 0;
  size_t j1 = 0, j2 = 0, j3 = 0;
  TLorentzVector vtau, vm12 ;
  double min_fvnC_2mu1tk = 10, min_fvnC = 100;
  int iTrk = 9999;
  vector<size_t> kinematic_muon_index;

  for(size_t i = 0; i < muons->size(); ++ i) 
    {
      const Muon & mu = (*muons)[i];
      if(!(mu.pt() > MuonPtCut_) || !(abs(mu.eta()) < MuonEtaCut_))continue;
      //    bool isID = false;
      if(mu.isPFMuon() && mu.isGlobalMuon()) kinematic_muon_index.push_back(i);// isID=true;
      //if(mu.isPFMuon() && (mu.isGlobalMuon()||mu.isTrackerMuon())) isID=true;
      
      //if( muon::isGoodMuon(mu, muon::TMOneStationTight)
      //  && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
      //  && mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0
      //  //if(!(m_1.innerTrack()->quality(TrackBase::highPurity)))continue;
      //  && abs(mu.innerTrack()->dxy(beamSpotHandle->position())) < 0.3
      //  && abs(mu.innerTrack()->dz(beamSpotHandle->position())) < 20
      //) isID=true;
      
      //    if(isID)kinematic_muon_index.push_back(i);
    }
  std::cout<<"number of muons  "<< kinematic_muon_index.size() << std::endl;
  if(kinematic_muon_index.size() < (do2mu_?2:3)) return;  //  if kinematic_muon < 3  - continue
  h_step->Fill(3);



  for(size_t i = 0; i < kinematic_muon_index.size()-1; ++ i)      //  loop over muons passed kinematic cuts
    {
    const Muon & m_1 = (*muons)[kinematic_muon_index[i]];


    for(size_t j = i+1; j < kinematic_muon_index.size(); ++ j)   // second loop to find a muon pair
      {
      const Muon & m_2 = (*muons)[kinematic_muon_index[j]];

      double dz12 = abs(m_2.innerTrack()->dz(beamSpotHandle->position())-m_1.innerTrack()->dz(beamSpotHandle->position()));  //   DeltaZ
      double dr12 = deltaR(m_1.eta(), m_1.phi(), m_2.eta(), m_2.phi());                                                      //   DeltaR

      if(dz12>0.5 ||  dr12>0.8)continue; // if Delta POCA_Z of two muon candiadate  > 0.5 cm or large deltaR  - skip the pair  candidate
      if(j<kinematic_muon_index.size()-1) 
	{
	for(size_t k = j+1; k < kinematic_muon_index.size(); ++ k) 
	  {
	  const Muon & m_3 = (*muons)[kinematic_muon_index[k]];

	  size_t n_muons_pt2p5 = 0;
	  if(m_1.pt()>2.5)n_muons_pt2p5++;
	  if(m_2.pt()>2.5)n_muons_pt2p5++;
	  if(m_3.pt()>2.5)n_muons_pt2p5++;
	  if(n_muons_pt2p5<2)continue;  //  ?? only to reduce the ntuple size ?

	  double dz23 = abs(m_3.innerTrack()->dz(beamSpotHandle->position())-m_2.innerTrack()->dz(beamSpotHandle->position()));
	  double dz31 = abs(m_3.innerTrack()->dz(beamSpotHandle->position())-m_1.innerTrack()->dz(beamSpotHandle->position()));
	  if(dz23>0.5 || dz31>0.5)continue;
	  double dr23 = deltaR(m_3.eta(), m_3.phi(), m_2.eta(), m_2.phi());
	  double dr31 = deltaR(m_3.eta(), m_3.phi(), m_1.eta(), m_1.phi());
	  if(dr23>0.8 || dr31>0.8)continue;


	  std::cout<<"muon1 charge "<<m_1.charge() << std::endl;
	  std::cout<<"muon2 charge "<<m_2.charge() << std::endl;
	  std::cout<<"muon3 charge "<<m_3.charge() << std::endl;
	  if(abs(m_1.charge()+m_2.charge()+m_3.charge())>1.1)continue;

	  vector<TransientTrack> t_trks;   
	  ESHandle<TransientTrackBuilder> theB;
	  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
	  TrackRef trk1 = m_1.innerTrack();
	  TrackRef trk2 = m_2.innerTrack();
	  TrackRef trk3 = m_3.innerTrack();
	  t_trks.push_back(theB->build(trk1));
	  t_trks.push_back(theB->build(trk2));
	  t_trks.push_back(theB->build(trk3));
	  KalmanVertexFitter kvf;
	  TransientVertex fv = kvf.vertex(t_trks);
	  std::cout<<"  if the vertex of three muons valid ? " << fv.isValid() << std::endl;
	  if(!fv.isValid()) continue; // eff ? 
	  double fvnC_tmp = fv.totalChiSquared()/fv.degreesOfFreedom();

	  //vtau.SetPxPyPzE(m_1.px()+m_2.px()+m_3.px(), m_1.py()+m_2.py()+m_3.py(), m_1.pz()+m_2.pz()+m_3.pz(), m_1.energy()+m_2.energy()+m_3.energy());
	  
	  if(n_reco==0)n_reco=3;
	  else {n_reco++;}
	  //cout<<"eveN: "<<eventN<<"  mass: "<<vtau.M()<<"  pt: "<<vtau.Pt()<<"  fvnC: "<<fvnC_tmp<<endl;
	  
	  //if(vtau.M() >  max_mtau){ // keep the max mass
	  if(fvnC_tmp < min_fvnC)  // check the loose quality of the SV
	    {
	      if(m_1.p()>m_2.p())  // sort muons by momentum
		{
		  if(m_2.p()>m_3.p())
		    {
		      j1=kinematic_muon_index[i]; j2=kinematic_muon_index[j]; j3=kinematic_muon_index[k];
		    }
		  else if(m_1.p()>m_3.p())
		    {
		      j1=kinematic_muon_index[i]; j2=kinematic_muon_index[k]; j3=kinematic_muon_index[j];
		    }
		  else
		    {
		      j1=kinematic_muon_index[k]; j2=kinematic_muon_index[i]; j3=kinematic_muon_index[j];
		    }
		}
	      else 
		{
		  if(m_1.p()>m_3.p())
		    {
		      j1=kinematic_muon_index[j]; j2=kinematic_muon_index[i]; j3=kinematic_muon_index[k];
		    }
		  else if(m_2.p()>m_3.p())
		    {
		      j1=kinematic_muon_index[j]; j2=kinematic_muon_index[k]; j3=kinematic_muon_index[i];
		    }
		  else 
		    {
		      j1=kinematic_muon_index[k]; j2=kinematic_muon_index[j]; j3=kinematic_muon_index[i];
		    }
		}
	      
	      min_fvnC = fvnC_tmp; // select finally the three_muon candidate by the best vertex;
	    }
	  }
	} 
      
      if(n_reco<3 && do2mu_ && m_1.pt() > 2.5 && m_2.pt() > 2.5)  ////////// if do 2mu+1trk
	{
	  for(size_t itk = 0; itk < trks->size(); itk++)
	  {
	    const Track & t = (*trks)[itk];

	    if(!isGoodTrack(t)) continue;
	    if(!(abs(t.dxy(beamSpotHandle->position())) < .3)  ||  !(abs(t.dz(beamSpotHandle->position())) < 20)) continue;   // check if the tracks is far from the BS.


	    
	    double dz23 = abs(t.dz(beamSpotHandle->position())-m_2.innerTrack()->dz(beamSpotHandle->position()));  // if the POCA of the track candidate is far from the muons - continue
	    double dz31 = abs(t.dz(beamSpotHandle->position())-m_1.innerTrack()->dz(beamSpotHandle->position()));

	    if(dz23 > 0.5 || dz31 > 0.5)  continue;

	    double dr23 = deltaR(t.eta(), t.phi(), m_2.eta(), m_2.phi());
	    double dr31 = deltaR(t.eta(), t.phi(), m_1.eta(), m_1.phi());

	    if(dr23 > 1.2 || dr31 > 1.2)    continue;
	    if(dr23 < 0.02 || dr31 < 0.02)  continue;

	    if( abs(m_1.charge()+m_2.charge()+t.charge())>1.1 ) continue;  // check the charge
	    
	    vector<TransientTrack> t_trks;
	    ESHandle<TransientTrackBuilder> theB;
	    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
	    TrackRef trk1 = m_1.innerTrack();
	    TrackRef trk2 = m_2.innerTrack();
	    TrackRef trk3 = TrackRef(trks, itk);
	    t_trks.push_back(theB->build(trk1));
	    t_trks.push_back(theB->build(trk2));
	    t_trks.push_back(theB->build(trk3));
	    KalmanVertexFitter kvf;
	    TransientVertex fv = kvf.vertex(t_trks);
	    if(!fv.isValid()) continue;
	    double fv_tC = fv.totalChiSquared();
	    double fv_dOF = fv.degreesOfFreedom();
	    double fv_nC = fv_tC/fv_dOF;
	    if(fv_nC > 5) continue;  // why 5 ?
	    if(fv_nC < min_fvnC_2mu1tk){
	      
	      //double t_energy = sqrt(0.140*0.140 + t.p()*t.p()); // pion mass
	      //vtau.SetPxPyPzE(m_1.px()+m_2.px()+t.px(), m_1.py()+m_2.py()+t.py(), m_1.pz()+m_2.pz()+t.pz(), m_1.energy()+m_2.energy()+t_energy);
	      //vm12.SetPxPyPzE(m_1.px()+m_2.px(), m_1.py()+m_2.py(), m_1.pz()+m_2.pz(), m_1.energy()+m_2.energy());
	      //if(min_fvnC_2mu1tk>9.99)cout<<endl;
	      //cout<<eventN<<"\t"<<fv_nC<<"\t"<<vm12.M()<<"\t"<<vtau.M()<<endl;
	      //if(vtau.M() >  max_mtau_2mu[0]t) // keep the max mass
	      
	      iTrk=itk; // index of the survived candidate
	      n_reco=2; //  2mu + track category
	      
	      if(m_1.p()>m_2.p()) // sort muons by p
		{
		  j1=kinematic_muon_index[i]; j2=kinematic_muon_index[j];
		}
	      else 
		{
		  j1=kinematic_muon_index[j]; j2=kinematic_muon_index[i];
		}

	      //max_mtau_2mu[0]t = vtau.M();
	      min_fvnC_2mu1tk = fv_nC;
	      
	    }
	    
	  } // loop of tracks
	
	} // if (n_reco<3)
      
      }
    
    }

  h_n3mu->Fill(n_reco);
  if(n_reco < (do2mu_?2:3)) return; 
  h_step->Fill(4);

  vector<Muon> mu; // -------------  this must be checked again
  mu.push_back((*muons)[j1]);
  mu.push_back((*muons)[j2]);
  mu.push_back((*muons)[j3]);


  ///////////////////
  // Trigger Objects
  edm::InputTag MuonFilterTag = edm::InputTag("hltTau3muTkVertexFilter", "", "HLT"); 
  size_t MuonFilterIndex = (*triggerSummary).filterIndex(MuonFilterTag); 
  trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();
  trigger::TriggerObjectCollection MuonLegObjects;
  if(MuonFilterIndex < (*triggerSummary).sizeFilters()) {
    const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(MuonFilterIndex);
    for(size_t j = 0; j < keysMuons.size(); j++ ){
      trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
      MuonLegObjects.push_back(foundObject);
    }
  }

  //cout<<MuonLegObjects.size()<<"\t"<<n_reco<<endl;

  if(n_reco >= 3 && doMC_){

    Handle<GenParticleCollection> genParticles2;
    iEvent.getByToken(genToken_, genParticles2);

    double dr = 0.03;

    //cout<<endl<<genParticles2->size()<<endl;
    for(size_t i = 2; i < genParticles2->size(); ++ i) {

      const GenParticle & p = (*genParticles2)[i];
      //const Candidate * mom = p.mother();

      //if(abs(p.pdgId())==13)cout<<p.status()<<" mom: "<<mom->pdgId()<<endl;
      if(p.charge()==0)continue;
      if(p.status()!=1)continue;
      if(p.p()<2.5)continue;
      if(abs(p.eta())>2.45)continue;
      //double dpt=3*(0.6+abs(p.eta()))/100;

      const Candidate * mom = p.mother();

      for(int j = 0; j < 3; j++) {

	if(pdgid_reco[j]==0){
	  double dpt = 3*mu[j].muonBestTrack()->ptError();
	  double dr1 = deltaR(mu[j].eta(), mu[j].phi(), p.eta(), p.phi());
	  double dpt1 = abs(mu[j].pt()-p.pt());
	  if(dr1 < dr && mu[j].charge()==p.charge() && dpt1 < dpt) {
	    pdgid_reco[j] = p.pdgId();
	    momid_reco[j] = mom->pdgId();
	    vxy_reco[j] = sqrt(p.vx()*p.vx() + p.vy()*p.vy());
	  }
	}

      }
    }
  }

  TrackRef t3 = TrackRef(trks, 0);
  if(n_reco>2) t3 = mu[2].innerTrack();
  else t3 = TrackRef(trks, iTrk);


  double t3_energy = sqrt(0.140*0.140 + t3->p()*t3->p()); // pion
  if(n_reco>2) t3_energy = sqrt(0.106*0.106 + t3->p()*t3->p()); // muon

  vtau.SetPxPyPzE(mu[0].px()+mu[1].px()+t3->px(), mu[0].py()+mu[1].py()+t3->py(), mu[0].pz()+mu[1].pz()+t3->pz(), mu[0].energy()+mu[1].energy()+t3_energy);

  m3mu_reco = vtau.M();

  double pt12 = (mu[0].pt()+mu[1].pt());
  double eta12 = (mu[0].eta()*mu[0].pt() + mu[1].eta()*mu[1].pt())/pt12;
  double phi12 = (mu[0].phi()*mu[0].pt() + mu[1].phi()*mu[1].pt())/pt12;
  m3mu_simp = sqrt(2*pt12*t3->pt()*( cosh(eta12-t3->eta()) - cos(phi12-t3->phi())));

  if(wideSB_) {
    if(m3mu_reco>3.2||m3mu_reco<1.)return;
  } else {
    if(m3mu_reco>2.05||m3mu_reco<1.6)return;
  }

  h_step->Fill(5);


  ////////////////////
  // 2 mu mass
  TLorentzVector mv1, mv2, mv3;
  vector<TLorentzVector> mv;
  mv1.SetPtEtaPhiM(mu[0].pt(), mu[0].eta(), mu[0].phi(), 0.106);
  mv2.SetPtEtaPhiM(mu[1].pt(), mu[1].eta(), mu[1].phi(), 0.106);
  mv3.SetPtEtaPhiM(t3->pt(), t3->eta(), t3->phi(), n_reco>2?0.106:0.140);
  m2mu_12 = (mv1+mv2).M();

  if(n_reco==2 && abs(m2mu_12-1.02)>0.02)return;

  m2mu_23 = (mv2+mv3).M();
  m2mu_31 = (mv3+mv1).M();
  m2mu_min= TMath::Min(m2mu_12, TMath::Min(m2mu_23, m2mu_31));
  m2mu_max= TMath::Max(m2mu_12, TMath::Max(m2mu_23, m2mu_31));

  if(abs(mu[0].charge()+mu[1].charge())>1) {
    mv.push_back(mv3); mv.push_back(mv1); mv.push_back(mv2);
  }
  else if (abs(mu[1].charge()+t3->charge())>1) {
    mv.push_back(mv1); mv.push_back(mv2); mv.push_back(mv3);
  }
  else if (abs(t3->charge()+mu[0].charge())>1) {
    mv.push_back(mv2); mv.push_back(mv3); mv.push_back(mv1);
  }
  double m2mu_osa=(mv[0]+mv[1]).M();
  double m2mu_osb=(mv[0]+mv[2]).M();
  m2mu_ss = (mv[1]+mv[2]).M();
  m2mu_os1 = (TMath::Max(m2mu_osa, m2mu_osb));
  m2mu_os2 = (TMath::Min(m2mu_osa, m2mu_osb));

  pt3mu_reco = vtau.Pt();
  p3mu_reco = vtau.P();
  eta3mu_reco = vtau.Eta();
  pt2mu_12 = sqrt((mu[0].px()+mu[1].px())*(mu[0].px()+mu[1].px()) + (mu[0].py()+mu[1].py())*(mu[0].py()+mu[1].py()));

  pt_max = 0; pt_min = 999;
  if(mu[0].pt()>pt_max)pt_max=mu[0].pt();
  if(mu[1].pt()>pt_max)pt_max=mu[1].pt();
  if(t3->pt()>pt_max)pt_max=t3->pt();
  if(mu[0].pt()<pt_min)pt_min=mu[0].pt();
  if(mu[1].pt()<pt_min)pt_min=mu[1].pt();
  if(t3->pt()<pt_min)pt_min=t3->pt();


  // Event based trigger matching
  trigmat_new = 0;
  for(size_t it = 0; it < MuonLegObjects.size(); it = it+3) {
    const trigger::TriggerObject & to1 = MuonLegObjects[it];
    const trigger::TriggerObject & to2 = MuonLegObjects[it+1];
    const trigger::TriggerObject & to3 = MuonLegObjects[it+2];

    bool mat1 = false, mat2 = false, mat3 = false;
    for(int i = 0; i < 2 ; i++) {
      if(deltaR(mu[i].eta(), mu[i].phi(), to1.eta(), to1.phi())<0.03 && abs(mu[i].pt()-to1.pt())/mu[i].pt()<0.1) mat1=true;
      if(deltaR(mu[i].eta(), mu[i].phi(), to2.eta(), to2.phi())<0.03 && abs(mu[i].pt()-to2.pt())/mu[i].pt()<0.1) mat2=true;
      if(n_reco>2)if(deltaR(mu[i].eta(), mu[i].phi(), to3.eta(), to3.phi())<0.03 && abs(mu[i].pt()-to3.pt())/mu[i].pt()<0.1) mat3=true;
    }

    if(n_reco>2) {
      if(deltaR(t3->eta(), t3->phi(), to1.eta(), to1.phi())<0.03 && abs(t3->pt()-to1.pt())/t3->pt()<0.1) mat1=true;
      if(deltaR(t3->eta(), t3->phi(), to2.eta(), to2.phi())<0.03 && abs(t3->pt()-to2.pt())/t3->pt()<0.1) mat2=true;
    }
    if(deltaR(t3->eta(), t3->phi(), to3.eta(), to3.phi())<0.03 && abs(t3->pt()-to3.pt())/t3->pt()<0.1) mat3=true;

    if(mat1 && mat2 && mat3) {trigmat_new = 1; break;}
  }

  //cout<<hlt_doublemu3_tau3mu<<"    "<<MuonLegObjects.size()<<endl<<endl;
  //cout<<endl<<lumiN<<endl;

  //cout<<"muon: "<<endl;  
  for(int i = 0; i < 3; i++) {  // loop muon

    //if(hlt_doublemu3_tau3mu)cout<<"reco mu"<<i<<"  "<<mu[i].pdgId()<<" "<<mu[i].pt()<<" "<<mu[i].eta()<<" "<<mu[i].phi()<<endl;

    // global muon variables
    p_out[i] = 0; eta_out[i] = -10;
    phi_out[i] = -10; p_glb[i] = 0;
    eta_glb[i] = -10; phi_glb[i] = -10;
    glbnC_reco[i] = 100; nOVMH_reco[i] = 0;
    mSWVH_reco[i] = 0, qprod_reco[i] = 0;
    outerchi2_reco[i] = 1000; vmuonhitcomb_reco[i] = 0;
    rpchits_reco[i] = 0;
    cschits_sta1[i] = 0;
    cscchi2_sta1[i] = 100;
    cschits_sta2[i] = 0;
    cscchi2_sta2[i] = 100;
    cscdxdz_sta1[i] = 0;
    cscdxdz_sta2[i] = 0;
    cscdydz_sta1[i] = 0;
    cscdydz_sta2[i] = 0;
    cscnsegm_sta1[i] = 0;
    cscnsegm_sta2[i] = 0;


    if(i==2 && n_reco<3) {

      uSta_reco[i] = 0;
      tKink_reco[i] = 0;
      gKink_reco[i] = 0;
      tRC2_reco[i] = 0;
      sRC2_reco[i] = 0;
      cLP_reco[i] = 0;
      cLM_reco[i] = 0;
      lDist_reco[i] = 0;
      gDEP_reco[i] = 0;
      tMat_reco[i] = 0;
      gTP_reco[i] = 0;
      calem_reco[i] = 0;
      calemS9_reco[i] = 0;
      calemS25_reco[i] = 0;
      calhad_reco[i] = 0;
      calhadS9_reco[i] = 0;
      nOMS_reco[i] = 99;
      nOM_reco[i] = 0;
      comp2d_reco[i] = 0;
      calocomp_reco[i] = 0;
      segmcomp_reco[i] = 0;
      segmcomp_0[i] = 0;
      segmcomp_1[i] = 0;
      segmcomp_2[i] = 0;
      segmcomp_3[i] = 0;
      //segmcomp_4[i] = 0;
      //segmcomp_5[i] = 0;
      //segmcomp_6[i] = 0;
      //segmcomp_7[i] = 0;
      pf_reco[i] = 0;
      rpcmu_reco[i] = 0;
      Iso_sumPt[i] = 0;
      Iso_nTr[i] = 0;
      Iso_emEt[i] = 0;
      Iso_hadEt[i] = 0;
      Iso_eVE[i] = 0;
      Iso_hVE[i] = 0;
      tma_reco[i] = 0;
      tmost_reco[i] = 0;
      tmosat_reco[i] = 0;
      tmlst_reco[i] = 0;
      tmlsat_reco[i] = 0;
      tmlsolpt_reco[i] = 0;
      tmlsoblpt_reco[i] = 0;
      timeatipinouterr_reco[i] = -1;

      pterr_reco[i] = t3->ptError()/t3->pt();
      trkhp_reco[i] = t3->quality(TrackBase::highPurity);
      charge_reco[i] = t3->charge();
      p_in[i] = t3->p();
      eta_in[i] = t3->eta();
      phi_in[i] = t3->phi();
      p_reco[i] = t3->p();
      pt_reco[i] = t3->pt();
      eta_reco[i] = t3->eta();
      phi_reco[i] = t3->phi();
      nOVPH_reco[i] = t3->hitPattern().numberOfValidPixelHits();
      iTvF_reco[i] = t3->validFraction();
      tLWM_reco[i] = t3->hitPattern().trackerLayersWithMeasurement();
      pLWM_reco[i] = t3->hitPattern().pixelLayersWithMeasurement();
      nOVTH_reco[i] = t3->hitPattern().numberOfValidTrackerHits();
      nOLTH_reco[i] = t3->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS);
      nOLTHin_reco[i] = t3->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
      nOLTHout_reco[i] = t3->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS);
      drtau_reco[i] = deltaR(t3->eta(), t3->phi(), vtau.Eta(), vtau.Phi());
      iTnC_reco[i] = t3->normalizedChi2();
      //	    cout<<" p_in  "  << p_in[i] << " i  "<<i << endl;
      trigmat_reco[i] = 0;
      for(size_t it = 0; it < MuonLegObjects.size(); it ++) {
	const trigger::TriggerObject & to = MuonLegObjects[it];
	if(deltaR(eta_reco[i], phi_reco[i], to.eta(), to.phi())<0.03 && abs(pt_reco[i]-to.pt())/pt_reco[i]<0.1)trigmat_reco[i] = to.id();  // was 0.05 and 0.3
      }

      continue;
    }


    charge_reco[i] = mu[i].charge();
    uSta_reco[i] = mu[i].combinedQuality().updatedSta;
    tKink_reco[i] = mu[i].combinedQuality().trkKink;
    gKink_reco[i] = TMath::Log(2+mu[i].combinedQuality().glbKink);
    tRC2_reco[i] = mu[i].combinedQuality().trkRelChi2;
    sRC2_reco[i] = mu[i].combinedQuality().staRelChi2;
    cLP_reco[i] = mu[i].combinedQuality().chi2LocalPosition;
    cLM_reco[i] = mu[i].combinedQuality().chi2LocalMomentum;
    lDist_reco[i] = mu[i].combinedQuality().localDistance;
    gDEP_reco[i] = mu[i].combinedQuality().globalDeltaEtaPhi;
    tMat_reco[i] = mu[i].combinedQuality().tightMatch;
    gTP_reco[i] = mu[i].combinedQuality().glbTrackProbability;
    calem_reco[i] = mu[i].calEnergy().em;
    calemS9_reco[i] = mu[i].calEnergy().emS9;
    calemS25_reco[i] = mu[i].calEnergy().emS25;
    calhad_reco[i] = mu[i].calEnergy().had;
    calhadS9_reco[i] = mu[i].calEnergy().hadS9;
    nOMS_reco[i] = mu[i].numberOfMatchedStations();
    nOM_reco[i] = mu[i].numberOfMatches(reco::Muon::SegmentArbitration);
    comp2d_reco[i] = muon::isGoodMuon(mu[i], muon::TM2DCompatibilityTight);
    calocomp_reco[i] = muon::caloCompatibility(mu[i]);
    segmcomp_reco[i] = muon::segmentCompatibility(mu[i]);
    //segmcomp_0[i] = scnew(mu[i], reco::Muon::SegmentAndTrackArbitration, 0);
    //segmcomp_1[i] = scnew(mu[i], reco::Muon::SegmentAndTrackArbitration, 1);
    //segmcomp_2[i] = scnew(mu[i], reco::Muon::SegmentAndTrackArbitration, 2);
    //segmcomp_3[i] = scnew(mu[i], reco::Muon::SegmentAndTrackArbitration, 3);
    //segmcomp_4[i] = scnew(mu[i], reco::Muon::SegmentAndTrackArbitration, 4);
    //segmcomp_5[i] = scnew(mu[i], reco::Muon::SegmentAndTrackArbitration, 5);
    //segmcomp_6[i] = scnew(mu[i], reco::Muon::SegmentAndTrackArbitration, 6);
    //segmcomp_7[i] = scnew(mu[i], reco::Muon::SegmentAndTrackArbitration, 7);

    //if(mu[i].pt()<8 && abs(mu[i].eta())<1.2 && mu[i].isGlobalMuon()) 
    //
    //if(mu[i].isGlobalMuon()) {
    //  cout<<endl<<"mu "<<i<<"  p: "<<mu[i].p()<<"  eta: "<<mu[i].eta()<<"  nOMS: "<<nOMS_reco[i]
    //      <<"  pdgid: "<<pdgid_reco[i]<<"  momid: "<<momid_reco[i]<<endl;
    //  cout<<"segmComp  Old: "<<segmcomp_reco[i]<<"  New: "<<scnew(mu[i], reco::Muon::SegmentAndTrackArbitration, 7)<<endl;
    //}
    trkhp_reco[i] = mu[i].innerTrack()->quality(TrackBase::highPurity);
    pf_reco[i] = mu[i].isPFMuon();
    rpcmu_reco[i] = mu[i].isRPCMuon();
    p_reco[i] = mu[i].p();
    pt_reco[i] = mu[i].pt();
    eta_reco[i] = mu[i].eta();
    phi_reco[i] = mu[i].phi();
    pterr_reco[i] = mu[i].muonBestTrack()->ptError()/mu[i].muonBestTrack()->pt();
    p_in[i] = mu[i].innerTrack()->p();
    eta_in[i] = mu[i].innerTrack()->eta();
    phi_in[i] = mu[i].innerTrack()->phi();
    nOVPH_reco[i] = mu[i].innerTrack()->hitPattern().numberOfValidPixelHits();
    iTvF_reco[i] = mu[i].innerTrack()->validFraction();
    tLWM_reco[i] = mu[i].innerTrack()->hitPattern().trackerLayersWithMeasurement();
    pLWM_reco[i] = mu[i].innerTrack()->hitPattern().pixelLayersWithMeasurement();
    nOVTH_reco[i] = mu[i].innerTrack()->hitPattern().numberOfValidTrackerHits();
    nOLTH_reco[i] = mu[i].innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS);
    nOLTHin_reco[i] = mu[i].innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
    nOLTHout_reco[i] = mu[i].innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS);
    drtau_reco[i] = deltaR(mu[i].eta(), mu[i].phi(), vtau.Eta(), vtau.Phi());
    Iso_sumPt[i] = mu[i].isolationR03().sumPt;
    Iso_nTr[i] = mu[i].isolationR03().nTracks;
    Iso_emEt[i] = mu[i].isolationR03().emEt;
    Iso_hadEt[i] = mu[i].isolationR03().hadEt;
    Iso_eVE[i] = mu[i].isolationR03().emVetoEt;
    Iso_hVE[i] = mu[i].isolationR03().hadVetoEt;
    iTnC_reco[i] = mu[i].innerTrack()->normalizedChi2();
    tma_reco[i] = muon::isGoodMuon(mu[i], muon::TrackerMuonArbitrated);
    tmost_reco[i] = muon::isGoodMuon(mu[i], muon::TMOneStationTight);
    tmosat_reco[i] = muon::isGoodMuon(mu[i], muon::TMOneStationAngTight);
    tmlst_reco[i] = muon::isGoodMuon(mu[i], muon::TMLastStationTight);
    tmlsat_reco[i] = muon::isGoodMuon(mu[i], muon::TMLastStationAngTight);
    tmlsolpt_reco[i] = muon::isGoodMuon(mu[i], muon::TMLastStationOptimizedLowPtTight);
    tmlsoblpt_reco[i] = muon::isGoodMuon(mu[i], muon::TMLastStationOptimizedBarrelLowPtTight);
    timeatipinouterr_reco[i] = mu[i].time().timeAtIpInOutErr;


    trigmat_reco[i] = 0;
    for(size_t it = 0; it < MuonLegObjects.size(); it ++) {
      const trigger::TriggerObject & to = MuonLegObjects[it];
      //if(i==0)cout<<"hlt  mu"<<it<<"  "<<to.id()<<" "<<to.pt()<<" "<<to.eta()<<" "<<to.phi()<<endl;
      if(deltaR(eta_reco[i], phi_reco[i], to.eta(), to.phi())<0.03 && abs(pt_reco[i]-to.pt())/pt_reco[i]<0.1)trigmat_reco[i] = to.id();  // was 0.05 and 0.3
    } 


    // global muon variables
    if(!mu[i].isGlobalMuon()) {
      continue;
    }

    p_out[i] = mu[i].outerTrack()->p();
    eta_out[i] = mu[i].outerTrack()->eta();
    phi_out[i] = mu[i].outerTrack()->phi();
    outerchi2_reco[i] = mu[i].outerTrack()->normalizedChi2();
    qprod_reco[i] = mu[i].outerTrack()->charge()*mu[i].innerTrack()->charge();

    p_glb[i] = mu[i].globalTrack()->p();
    eta_glb[i] = mu[i].globalTrack()->eta();
    phi_glb[i] = mu[i].globalTrack()->phi();
    glbnC_reco[i] = mu[i].globalTrack()->normalizedChi2();
    nOVMH_reco[i] = mu[i].globalTrack()->hitPattern().numberOfValidMuonHits();
    mSWVH_reco[i] = mu[i].outerTrack()->hitPattern().muonStationsWithValidHits();

    unsigned int dt1(0),dt2(0),dt3(0),dt4(0);
    unsigned int rpc1(0),rpc2(0),rpc3(0),rpc4(0);
    unsigned int csc1(0),csc2(0),csc3(0),csc4(0);
    double comb(0);
    const reco::HitPattern &pattern = mu[i].outerTrack()->hitPattern();
    for (int i=0;i<pattern.numberOfAllHits(reco::HitPattern::TRACK_HITS);i++)
      { 
	uint32_t hit = pattern.getHitPattern(reco::HitPattern::TRACK_HITS,i);
	if (pattern.validHitFilter(hit) != 1) {continue;}
	if (pattern.getMuonStation(hit) == 1)
	  { 
	    if (pattern.muonDTHitFilter(hit)) {dt1++;}
	    if (pattern.muonRPCHitFilter(hit)) {rpc1++;}
	    if (pattern.muonCSCHitFilter(hit)) {csc1++;}
	  }
	else if (pattern.getMuonStation(hit) == 2)
	  { 
	    if (pattern.muonDTHitFilter(hit)) {dt2++;}
	    if (pattern.muonRPCHitFilter(hit)) {rpc2++;}
	    if (pattern.muonCSCHitFilter(hit)) {csc2++;}
	  }
	else if (pattern.getMuonStation(hit) == 3)
	  { 
	    if (pattern.muonDTHitFilter(hit)) {dt3++;}
	    if (pattern.muonRPCHitFilter(hit)) {rpc3++;}
	    if (pattern.muonCSCHitFilter(hit)) {csc3++;}
	  }
	else if (pattern.getMuonStation(hit) == 4)
	  { 
	    if (pattern.muonDTHitFilter(hit)) {dt4++;}
	    if (pattern.muonRPCHitFilter(hit)) {rpc4++;}
	    if (pattern.muonCSCHitFilter(hit)) {csc4++;}
	  }    
      }
    comb = (dt1+dt2+dt3+dt4)/2. + (rpc1+rpc2+rpc3+rpc4);
    csc1>6 ? comb+=6 : comb+=csc1;
    csc2>6 ? comb+=6 : comb+=csc2;
    csc3>6 ? comb+=6 : comb+=csc3;
    csc4>6 ? comb+=6 : comb+=csc4;
    vmuonhitcomb_reco[i] = comb;
    rpchits_reco[i] = rpc1+rpc2+rpc3+rpc4;

  } // loop muon


  // some handy max, min variables
  drtau_max = TMath::Max(drtau_reco[0], TMath::Max(drtau_reco[1], drtau_reco[2] ));
  tLWM_min = TMath::Min(tLWM_reco[0], TMath::Min(tLWM_reco[1], tLWM_reco[2]));
  nOMS_min = TMath::Min(nOMS_reco[0], TMath::Min(nOMS_reco[1], nOMS_reco[2]));
  eta_min = TMath::Min(abs(eta_reco[0]), TMath::Min(abs(eta_reco[1]), abs(eta_reco[2])));
  eta_max = TMath::Max(abs(eta_reco[0]), TMath::Max(abs(eta_reco[1]), abs(eta_reco[2])));

  //MuonShower muonShowerInformation1 = (*muonShowerInformationValueMapH_)[muonRef1];
  //nSCH_reco[0] = muonShowerInformation1.nStationCorrelatedHits.at(1);

  dr12_reco = deltaR(mu[0].eta(), mu[0].phi(), mu[1].eta(), mu[1].phi());
  dr23_reco = deltaR(t3->eta(), t3->phi(), mu[1].eta(), mu[1].phi());
  dr31_reco = deltaR(mu[0].eta(), mu[0].phi(), t3->eta(), t3->phi());

  dr_min = 99;
  if(dr12_reco<dr_min) {dr_min = dr12_reco; ifar=2;}
  if(dr23_reco<dr_min) {dr_min = dr23_reco; ifar=0;}
  if(dr31_reco<dr_min) {dr_min = dr31_reco; ifar=1;}


  ////////////////////
  // Fit 3mu vertex
  vector<TransientTrack> t_trks;
  TrackRef trk1 = mu[0].innerTrack();
  TrackRef trk2 = mu[1].innerTrack();
  ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  t_trks.push_back(theB->build(trk1));
  t_trks.push_back(theB->build(trk2));
  t_trks.push_back(theB->build(t3));


  //TrajectoryStateClosestToPoint mu1TS = t_trks[0].impactPointTSCP();
  //TrajectoryStateClosestToPoint mu2TS = t_trks[1].impactPointTSCP();
  ClosestApproachInRPhi cApp12, cApp23, cApp31;
  cApp12.calculate(t_trks[0].initialFreeState(), t_trks[1].initialFreeState());
  cApp23.calculate(t_trks[1].initialFreeState(), t_trks[2].initialFreeState());
  cApp31.calculate(t_trks[2].initialFreeState(), t_trks[0].initialFreeState());
  //cApp.calculate(mu1TS.theState(), mu2TS.theState());
  if(!(cApp12.status()&&cApp23.status()&&cApp31.status())) { return; cout<<"DCA unvalid!"<<endl; }
  dca12_reco = cApp12.distance();
  dca23_reco = cApp23.distance();
  dca31_reco = cApp31.distance();
  dca_max = TMath::Max(dca12_reco, TMath::Max(dca23_reco, dca31_reco));

  KalmanVertexFitter kvf(true);
  TransientVertex fv = kvf.vertex(t_trks);
  if(!fv.isValid()) { return; cout<<"Vertex Fit unvalid!"<<endl; }

  TLorentzVector vtau_refit, vmu_refit;
  vtau_refit.SetPtEtaPhiM(0, 0, 0, 0);
  vector<TransientTrack>::const_iterator trkIt = fv.refittedTracks().begin();
  for(; trkIt != fv.refittedTracks().end(); ++ trkIt) {
    const Track & trkrefit = trkIt->track();
    vmu_refit.SetPtEtaPhiM(trkrefit.pt(), trkrefit.eta(), trkrefit.phi(), (n_reco>2?0.106:0.140));
    vtau_refit += vmu_refit;
  }
  m3mu_refit = vtau_refit.M();

  fv_tC = fv.totalChiSquared();
  int fv_dOF = fv.degreesOfFreedom();
  fv_nC = fv_tC/fv_dOF;
  fv_Prob = TMath::Prob(fv_tC,(int)fv_dOF);

  //GlobalError err = fv.positionError(); //(verr.At(0,0), verr.At(1,0), verr.At(1,1), verr.At(2,0), verr.At(2,1), verr.At(2,2) );
  //GlobalPoint displacementFromBeamspot( -1*((bs.x0() -  fv.position().x()) +  (fv.position().z() - bs.z0()) * bs.dxdz()),
  //-1*((bs.y0() - fv.position().y())+  (fv.position().z() - bs.z0()) * bs.dydz()), 0);
  //fv_dxy = displacementFromBeamspot.perp();
  //double fv_dxyerr = sqrt(err.rerr(displacementFromBeamspot));
  //fv_dxysig = fv_dxy/fv_dxyerr;

  vector<TransientTrack> t_trks12, t_trks23, t_trks31;
  t_trks12.push_back(theB->build(trk1)); t_trks12.push_back(theB->build(trk2));
  t_trks23.push_back(theB->build(trk2)); t_trks23.push_back(theB->build(t3));
  t_trks31.push_back(theB->build(t3)); t_trks31.push_back(theB->build(trk1));
  KalmanVertexFitter kvf_trks12, kvf_trks23, kvf_trks31;
  TransientVertex fv_trks12 = kvf_trks12.vertex(t_trks12);
  TransientVertex fv_trks23 = kvf_trks23.vertex(t_trks23);
  TransientVertex fv_trks31 = kvf_trks31.vertex(t_trks31);
  fvwo_tC[0] = fv_trks23.totalChiSquared();
  fvwo_nC[0] = fvwo_tC[0]/fv_trks23.degreesOfFreedom();
  fvwo_tC[1] = fv_trks31.totalChiSquared();
  fvwo_nC[1] = fvwo_tC[1]/fv_trks31.degreesOfFreedom();
  fvwo_tC[2] = fv_trks12.totalChiSquared();
  fvwo_nC[2] = fvwo_tC[2]/fv_trks12.degreesOfFreedom();


  TVector3 vtauxyz(vtau.Px(), vtau.Py(), vtau.Pz());

  ////////////////////
  // find the good PV
  double PVZ = fv.position().z()-fv_dxy*vtau.Pz()/vtau.Pt();
  double dispv1 = 99, dispvgen=99, dphi_pv = -1;
  //int ipvPVZ = 99, ipvgen = 99;
  for(size_t jpv = 0; jpv < pvs->size(); jpv++) {
    const Vertex & vi = (*pvs)[jpv];

    if(abs(vi.position().Z()-PVZ)<dispv1){
      dispv1=abs(vi.position().Z()-PVZ);
      ipv1=jpv;
    }
    if(abs(vi.position().Z()-gen_pv)<dispvgen){
      dispvgen=abs(vi.position().Z()-gen_pv);
      ipv_gen=jpv;
    }
    TVector3 Dv3D_reco(fv.position().x() - vi.x(), fv.position().y() - vi.y(), fv.position().z() - vi.z());
    double Cosdphi_3D = Dv3D_reco.Dot(vtauxyz)/(Dv3D_reco.Mag()*vtauxyz.Mag());
    if(Cosdphi_3D>dphi_pv){
      dphi_pv = Cosdphi_3D;
      ipv2=jpv;
    }

  }
  const Vertex & pv0 = (*pvs)[ipv2];


  //////////////////////////////////////////////////
  // refit PV with and w.o. the 3 mu
  pv1_tC = 999; pv1_nC = 99; pv2_tC = 999; pv2_nC = 99;
  vector<TransientTrack> pv_trks;
  TransientVertex pv2, pv1;

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  for(Vertex::trackRef_iterator itk = pv0.tracks_begin(); itk != pv0.tracks_end(); itk++) {
    if((**itk).pt()>1) {
      if(deltaR(mu[0].eta(), mu[0].phi(), (**itk).eta(), (**itk).phi())<0.01)continue;
      if(deltaR(mu[1].eta(), mu[1].phi(), (**itk).eta(), (**itk).phi())<0.01)continue;
      if(deltaR(t3->eta(), t3->phi(), (**itk).eta(), (**itk).phi())<0.01)continue;
    }
    pv_trks.push_back(theB->build(**itk));
  }
  if(pv_trks.size()>1) {
    KalmanVertexFitter kvf_pv;
    pv1 = kvf_pv.vertex(pv_trks);
    if(pv1.isValid()){
      pv1_tC = pv1.totalChiSquared();
      pv1_nC = pv1_tC/pv1.degreesOfFreedom();
    }

    // adding the 3 mu tracks
    pv_trks.push_back(theB->build(trk1));
    pv_trks.push_back(theB->build(trk2));
    pv_trks.push_back(theB->build(t3));
    pv2 = kvf_pv.vertex(pv_trks);
    if(pv2.isValid()){
      pv2_tC = pv2.totalChiSquared();
      pv2_nC = pv2_tC/pv2.degreesOfFreedom();
    }
  }
  //cout<<pv0.x()-pv1.position().x()<<"  "<<pv0.y()-pv1.position().y()<<"  "<<pv0.z()-pv1.position().z()<<"  "<<endl;

  Vertex pvv = pv0;  // the final PV
  if(pv1.isValid()) pvv = Vertex(pv1);
  math::XYZPoint pv1P = math::XYZPoint(pvv.x(), pvv.y(), pvv.z());


  d0_reco[0] = abs(mu[0].innerTrack()->dxy(pv1P));
  d0_reco[1] = abs(mu[1].innerTrack()->dxy(pv1P));
  d0_reco[2] = abs(t3->dxy(pv1P));
  d0sig_reco[0] = -1; d0sig_reco[1] = -1; d0sig_reco[2] = -1;
  GlobalVector dir1(mu[0].px(), mu[0].py(), mu[0].pz());
  GlobalVector dir2(mu[1].px(), mu[1].py(), mu[1].pz());
  GlobalVector dir3(t3->px(), t3->py(), t3->pz());
  std::pair<bool, Measurement1D> ip2d_1 = IPTools::signedTransverseImpactParameter(t_trks[0], dir1, pvv);
  std::pair<bool, Measurement1D> ip2d_2 = IPTools::signedTransverseImpactParameter(t_trks[1], dir2, pvv);
  std::pair<bool, Measurement1D> ip2d_3 = IPTools::signedTransverseImpactParameter(t_trks[2], dir3, pvv);
  if(ip2d_1.first) d0sig_reco[0] = abs(ip2d_1.second.value()/ip2d_1.second.error());
  if(ip2d_2.first) d0sig_reco[1] = abs(ip2d_2.second.value()/ip2d_2.second.error());
  if(ip2d_3.first) d0sig_reco[2] = abs(ip2d_3.second.value()/ip2d_3.second.error());


  ////////////////////
  // displacement 2D
  TVector3 dv2D_reco(fv.position().x() - pv1P.x(), fv.position().y() - pv1P.y(), 0);
  TVector3 vtauxy(vtau.Px(), vtau.Py(), 0);
  fv_cosdphi = dv2D_reco.Dot(vtauxy)/(dv2D_reco.Perp()*vtauxy.Perp());
  VertexDistanceXY vdistXY;
  Measurement1D distXY = vdistXY.distance(Vertex(fv), pvv);
  fv_dxy = distXY.value();
  fv_dxysig = distXY.significance();
  fv_ppdl = distXY.value()*fv_cosdphi * m3mu_reco/vtauxy.Perp();


  ////////////////////
  // displacement 3D
  TVector3 dv3D_reco(fv.position().x() - pv1P.x(), fv.position().y() - pv1P.y(), fv.position().z() - pv1P.z());
  //TVector3 vtauxyz(vtau.Px(), vtau.Py(), vtau.Pz());
  fv_cosdphi3D = dv3D_reco.Dot(vtauxyz)/(dv3D_reco.Mag()*vtauxyz.Mag());
  VertexDistance3D dist;
  fv_d3D = dist.distance(Vertex(fv), pvv).value(); // = dv_reco.Mag() ??
  fv_d3Dsig = dist.distance(Vertex(fv), pvv).significance();
  fv_ppdl3D = fv_d3D*fv_cosdphi3D*m3mu_reco/vtau.P();



  vector<double> softmueta, softmuphi;
  pv_nmu = 0;

  for(size_t i = 0; i < muons->size(); ++ i) {
    if(i==j1)continue;
    if(i==j2)continue;
    if(i==j3)continue;
    const Muon & m_1 = (*muons)[i];
    if(!(abs(m_1.eta())<2.4)) continue;
    if(!(muon::isGoodMuon(m_1, muon::TMOneStationTight))) continue;
    if(!(m_1.innerTrack()->hitPattern().trackerLayersWithMeasurement()>5))continue;
    if(!(m_1.innerTrack()->hitPattern().pixelLayersWithMeasurement()>0))continue;
    //if(!(abs(m_1.innerTrack()->dxy(pv0.position())) < .3))continue;
    if(!(abs(m_1.innerTrack()->dz(pv1P)) < 1))continue;
    pv_nmu++;
    softmueta.push_back(m_1.eta());
    softmuphi.push_back(m_1.phi());
  }

  //for(Vertex::trackRef_iterator itk = pv0.tracks_begin(); itk != pv0.tracks_end(); itk++) {
  //  for(size_t imu = 0; imu < softmueta.size(); imu ++) {
  //    if(deltaR(softmueta[imu], softmuphi[imu], (**itk).eta(), (**itk).phi())<0.01)
  //      pv_nmu++;
  //  }
  //}

  ////////////////////
  // secondary vertices
  n_sv = 0 ;
  for(size_t isv = 0; isv < svs->size(); isv++) {
    const Vertex & sv = (*svs)[isv];
    if(abs(sv.p4().M()-0.498)<.03 && sv.tracksSize()==2)continue; // no Ks

    double dx = sv.x()-pv1P.x();
    double dy = sv.y()-pv1P.y();
    double dz = sv.z()-pv1P.z();
    TVector3 sv_reco(dx, dy, dz);
    sv_overlap[n_sv]=deltaR(sv_reco.Eta(), sv_reco.Phi(), dv3D_reco.Eta(), dv3D_reco.Phi());

    TVector3 svxyz(sv.p4().Px(), sv.p4().Py(), sv.p4().Pz());
    sv_cosdphi3D[n_sv] = sv_reco.Dot(svxyz)/(sv_reco.Mag()*svxyz.Mag());
    VertexDistance3D distsv;
    sv_d3D[n_sv] = distsv.distance(sv, pvv).value();
    sv_d3Dsig[n_sv] = distsv.distance(sv, pvv).significance();
    sv_ppdl3D[n_sv] = sv_d3D[n_sv]*sv_cosdphi3D[n_sv]*sv.p4().M()/sv.p4().P();

    sv_nmu[n_sv] = 0;
    for(Vertex::trackRef_iterator itk = sv.tracks_begin(); itk != sv.tracks_end(); itk++) {
      for(size_t imu = 0; imu < softmueta.size(); imu ++) {
	if(deltaR(softmueta[imu], softmuphi[imu], (**itk).eta(), (**itk).phi())<0.01)
	  sv_nmu[n_sv] ++;
      }
    }

    sv_mass[n_sv] = sv.p4().M();
    sv_pt[n_sv] = sv.p4().Pt();
    sv_dz[n_sv] = abs(dz);
    sv_ntrk[n_sv] = sv.tracksSize();
    n_sv++;
  }

  dzpv_reco[0] = 1;
  dzpv_reco[1] = 1;
  dzpv_reco[2] = 1;

  dz_reco[0] = abs(mu[0].innerTrack()->dz(pv1P));
  dz_reco[1] = abs(mu[1].innerTrack()->dz(pv1P));
  dz_reco[2] = abs(t3->dz(pv1P));

  for(size_t jpv = 0; jpv < pvs->size(); jpv++) {
    if(jpv==ipv2)continue;
    const Vertex & vi = (*pvs)[jpv];
    if(abs(mu[0].innerTrack()->dz(vi.position()))<dz_reco[0]) dzpv_reco[0]=-1;
    if(abs(mu[1].innerTrack()->dz(vi.position()))<dz_reco[1]) dzpv_reco[1]=-1;
    if(abs(t3->dz(vi.position()))<dz_reco[2]) dzpv_reco[2]=-1;
  }


  ////////////////////
  // Track Isolation
  // How to decide if a track is associated with a certain PV ?
  double pttrk_tau = 0, pttrk_tau05 = 0,  pttrk_m1 = 0, pttrk_m2 = 0, pttrk_m3 = 0;
  mindca_iso = 99; mindca_iso05 = 99;
  ntrk_tau = 0; ntrk_tau05 = 0; ntrk_tau_b = 0; ntrk_sum = 0;
  ntrk_reco[0] = 0;  ntrk_reco[1] = 0;  ntrk_reco[2] = 0;
  ntrk0p1 = 0; ntrk0p2 = 0; ntrk0p5 = 0; maxdxy_pv0 =0;

  math::XYZPoint fvP = math::XYZPoint(fv.position().x(), fv.position().y(), fv.position().z());
  for(size_t i = 0; i < trks->size(); i++) {
    const Track & t = (*trks)[i];
    if(!(t.quality(TrackBase::tight)))continue;
    if(deltaR(mu[0].eta(), mu[0].phi(), t.eta(), t.phi())<0.01)continue;
    if(deltaR(mu[1].eta(), mu[1].phi(), t.eta(), t.phi())<0.01)continue;
    if(deltaR(t3->eta(), t3->phi(), t.eta(), t.phi())<0.01)continue;

    double dz = abs(t.dz(fvP));
    double dxy = abs(t.dxy(fvP));
    double dca_fv = sqrt(dz*dz+dxy*dxy);
    double dr_tau = deltaR(t.eta(), t.phi(), vtau.Eta(), vtau.Phi());

    // iso no. 1b - using pt_min, drtau_max of the 3 mu
    if(t.pt() > 0.33*pt_min && dr_tau < 3.*drtau_max && dca_fv<0.05 ) {
      pttrk_tau += t.pt();
      ntrk_tau++; // iso 3b
      if(dca_fv<mindca_iso)mindca_iso=dca_fv; // iso 4b
    } 

    if(t.pt()<1.0) continue;  // was 1.2
    // iso no. 1
    if(dr_tau < 0.5 && dca_fv<0.05 ) {
      pttrk_tau05 += t.pt();
      ntrk_tau05++; // iso 3
      //if(dca_fv<mindca_iso05)mindca_iso05=dca_fv; // iso 4
    }

    if(dca_fv<0.05)ntrk_tau_b++; // iso 3b
    if(dca_fv<mindca_iso05)mindca_iso05=dca_fv; // iso 4


    TransientTrack trkiso = theB->build(t);
    ClosestApproachInRPhi cAppm1, cAppm2, cAppm3;
    cAppm1.calculate(trkiso.initialFreeState(), t_trks[0].initialFreeState());
    cAppm2.calculate(trkiso.initialFreeState(), t_trks[1].initialFreeState());
    cAppm3.calculate(trkiso.initialFreeState(), t_trks[2].initialFreeState());
    if(!(cAppm1.status()&&cAppm2.status()&&cAppm3.status())) continue;

    // iso no. 2
    if(deltaR(t.eta(), t.phi(), mu[0].eta(), mu[0].phi()) < 0.3 && cAppm1.distance() < 0.1) {// && dz1 < .3) 
      ntrk_reco[0]++;
      pttrk_m1 += t.pt();
    }
    if(deltaR(t.eta(), t.phi(), mu[1].eta(), mu[1].phi()) < 0.3 && cAppm2.distance() < 0.1) {//&& dz2 < .3) 
      ntrk_reco[1]++;
      pttrk_m2 += t.pt();
    }
    if(deltaR(t.eta(), t.phi(), t3->eta(), t3->phi()) < 0.3 && cAppm3.distance() < 0.1) {//&& dz3 < .3) 
      ntrk_reco[2]++;
      pttrk_m3 += t.pt();
    }
    if( (deltaR(t.eta(), t.phi(), mu[0].eta(), mu[0].phi()) < 0.3 && cAppm1.distance() < 0.1 )
	||(deltaR(t.eta(), t.phi(), mu[1].eta(), mu[1].phi()) < 0.3 && cAppm2.distance() < 0.1 )
	||(deltaR(t.eta(), t.phi(), t3->eta(), t3->phi()) < 0.3 && cAppm3.distance() < 0.1 )
	) ntrk_sum++;


    // displaced track counting
    // only tracks consistent with PV
    double dz_pv0=abs(t.dz(pv1P));
    if(!(dz_pv0 < 1))continue;
    double dxy_pv0 = abs(t.dxy(pv1P));
    if(dxy_pv0>0.1) ntrk0p1++;
    if(dxy_pv0>0.2) ntrk0p2++;
    if(dxy_pv0>0.5) ntrk0p5++;
    if(dxy_pv0>maxdxy_pv0) maxdxy_pv0 = dxy_pv0;

  }

  trkrel_tau = pttrk_tau/vtau.Pt();
  trkrel_tau05 = pttrk_tau05/vtau.Pt();
  trkrel_reco[0] = pttrk_m1/mu[0].pt(); trkrel_reco[1] = pttrk_m2/mu[1].pt(); trkrel_reco[2] = pttrk_m3/t3->pt();
  trkrel_max = TMath::Max(trkrel_reco[0], TMath::Max( trkrel_reco[1], trkrel_reco[2]));


  // Good Global Muon, Tight Global Muon
  ggm_reco[0] = (glbnC_reco[0]<3 && cLP_reco[0]<12 && tKink_reco[0]<20 && segmcomp_reco[0]>0.303)?1:0;
  tgm_reco[0] = (glbnC_reco[0]<10 && pf_reco[0] && nOVMH_reco[0]>0 && nOMS_reco[0]>1
		 && d0_reco[0]<0.2 && dz_reco[0]<0.5 && nOVPH_reco[0]>0 && tLWM_reco[0]>5)?1:0;
  ggm_reco[1] = (glbnC_reco[1]<3 && cLP_reco[1]<12 && tKink_reco[1]<20 && segmcomp_reco[1]>0.303)?1:0;
  tgm_reco[1] = (glbnC_reco[1]<10 && pf_reco[1] && nOVMH_reco[1]>0 && nOMS_reco[1]>1
		 && d0_reco[1]<0.2 && dz_reco[1]<0.5 && nOVPH_reco[1]>0 && tLWM_reco[1]>5)?1:0;
  ggm_reco[2] = (glbnC_reco[2]<3 && cLP_reco[2]<12 && tKink_reco[2]<20 && segmcomp_reco[2]>0.303)?1:0;
  tgm_reco[2] = (glbnC_reco[2]<10 && pf_reco[2] && nOVMH_reco[2]>0 && nOMS_reco[2]>1
		 && d0_reco[2]<0.2 && dz_reco[2]<0.5 && nOVPH_reco[2]>0 && tLWM_reco[2]>5)?1:0;




  /////////////////
  // b tag
  //  to add Jet ID ??
  njet20 = 0;
  for(size_t j = 0 ; j < btagsCvsB->size(); j++) {
    const JetTag & btag1 = (*btagsCvsB)[j];

    if(btag1.first->pt()<20) break;
    if(deltaR(vtau.Eta(), vtau.Phi(), btag1.first->eta(), btag1.first->phi())<0.4)jet_overlap[njet20]=1;
    else jet_overlap[njet20]=0;

    jet_pt[njet20] = btag1.first->pt();
    btagcvsb[njet20] = btag1.second;
    //cout<<"jet"<<j<<"  pt: "<<btag.first->pt()<<"  btag: "<<btag.second<<endl;
    const JetTag & btag2 = (*btagsMVA)[j];
    btagmva[njet20] = btag2.second;
    const JetTag & btag3 = (*btagsCSV)[j];
    btagcsv[njet20] = (btag3.second<0 ? 0:btag3.second);

    njet20 ++;
  }

  tree->Fill();
  h_step->Fill(6);

}


void T3MNtuple::fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Handle<TrackCollection> trackCollection;
  iEvent.getByToken(trackToken_, trackCollection);

  std::vector<reco::Track>::const_iterator trIt  = trackCollection->begin();
  std::vector<reco::Track>::const_iterator trEnd = trackCollection->end();

  for (; trIt != trEnd; ++trIt) 
    {
      std::vector<double> iTrack_p4;
      std::vector<double> iTrack_poca;
      const reco::Track tr = (*trIt);
      iTrack_p4.push_back(sqrt(pow(tr.p(),2.0) + pow(PDGInfo::pi_mass(),2.0)));
      iTrack_p4.push_back(tr.px());
      iTrack_p4.push_back(tr.py());
      iTrack_p4.push_back(tr.pz());
      Track_p4.push_back(iTrack_p4);

      Track_normalizedChi2.push_back(tr.normalizedChi2());
      Track_numberOfValidHits.push_back(tr.numberOfValidHits());
      Track_charge.push_back(tr.charge());
      Track_dxy.push_back(tr.dxy());
      Track_dz.push_back(tr.dz());
      iTrack_poca.push_back(tr.vx());
      iTrack_poca.push_back(tr.vy());
      iTrack_poca.push_back(tr.vz());
      Track_poca.push_back(iTrack_poca);

      Track_dxyError.push_back(tr.dxyError());
      Track_dzError.push_back(tr.dzError());
    }
}

void T3MNtuple::fillMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Handle<MuonCollection> muonCollection;
  iEvent.getByToken(muonToken_, muonCollection);
  int Muon_index = 0;
  for (reco::MuonCollection::const_iterator iMuon = muonCollection->begin(); iMuon != muonCollection->end(); ++iMuon, Muon_index++) {
    reco::MuonRef RefMuon(muonCollection, Muon_index);
    //    if (isGoodMuon(RefMuon)) {
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


      Muon_numberOfChambers.push_back(RefMuon->numberOfChambers());
      Muon_isGlobalMuon.push_back(RefMuon->isGlobalMuon());
      Muon_isPFMuon.push_back(RefMuon->isPFMuon());
      Muon_isStandAloneMuon.push_back(RefMuon->isStandAloneMuon());
      Muon_isTrackerMuon.push_back(RefMuon->isTrackerMuon());
      Muon_isCaloMuon.push_back(RefMuon->isCaloMuon());
      Muon_isQualityValid.push_back(RefMuon->isQualityValid());
      Muon_isTimeValid.push_back(RefMuon->isTimeValid());
      Muon_isIsolationValid.push_back(RefMuon->isIsolationValid());
      Muon_numberOfMatchedStations.push_back(RefMuon->numberOfMatchedStations());
      Muon_numberOfMatches.push_back(RefMuon->numberOfMatches());
      Muon_charge.push_back(RefMuon->charge());


      if (RefMuon->isGlobalMuon()) {
	Muon_normChi2.push_back(RefMuon->globalTrack()->normalizedChi2());
	Muon_hitPattern_numberOfValidMuonHits.push_back(RefMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
	Muon_trackerLayersWithMeasurement.push_back(RefMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
	Muon_numberofValidPixelHits.push_back(RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits());

      } else {
	Muon_normChi2.push_back(0);
	Muon_hitPattern_numberOfValidMuonHits.push_back(0);
	Muon_trackerLayersWithMeasurement.push_back(0);
	Muon_numberofValidPixelHits.push_back(0);
      }
      if (RefMuon->isTrackerMuon()) {
	Muon_innerTrack_numberofValidHits.push_back(RefMuon->innerTrack()->numberOfValidHits());
	Muon_hitPattern_pixelLayerwithMeas.push_back(RefMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());

      } else {
	Muon_innerTrack_numberofValidHits.push_back(0);
	Muon_hitPattern_pixelLayerwithMeas.push_back(0);
      }



  }
}

void T3MNtuple::fillMCTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
if (!iEvent.isRealData())
  {
  Handle<vector<PileupSummaryInfo> >  PupInfo;
  iEvent.getByToken(puToken_, PupInfo);
  puN = PupInfo->begin()->getTrueNumInteractions();
  }
}


void T3MNtuple::fillL1(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  gtUtil_->retrieveL1(iEvent, iSetup, algToken_);
  const vector<pair<string, bool> > initialDecisions = gtUtil_->decisionsInitial();

  if (!iEvent.isRealData())
    {
      //      gtUtil_->retrieveL1(iEvent, iSetup, algToken_);
      //      const vector<pair<string, bool> > initialDecisions = gtUtil_->decisionsInitial();
      for (size_t i_l1t = 0; i_l1t < initialDecisions.size(); i_l1t++) 
	{
    	  string l1tName = (initialDecisions.at(i_l1t)).first;
	  if( l1tName == "NULL") continue;
    
	  if( l1tName == "L1_DoubleMu0" ) 
	    {
	      if( initialDecisions.at(i_l1t).second ) l1_doublemu0 = 1;
	    }
	  if( l1tName == "L1_TripleMu0" )
	    {
	      if( initialDecisions.at(i_l1t).second ) l1_triplemu0 = 1;
	    }
	  if( l1tName == "L1_TripleMu_5_0_0" )
	    {
	      if( initialDecisions.at(i_l1t).second ) l1_triplemu500 = 1;
	    }
	  if( l1tName == "L1_DoubleMu0er1p6_dEta_Max1p8" )
	    {
	      if( initialDecisions.at(i_l1t).second ) l1_doublemu0_eta1p6 = 1;
	    }
	  if( l1tName == "L1_DoubleMu0er1p6_dEta_Max1p8_OS" )
	    {
	      if( initialDecisions.at(i_l1t).second ) l1_doublemu0_eta1p6_os = 1;
	    }
	  if( l1tName == "L1_DoubleMu0er1p4_dEta_Max1p8_OS" ) 
	    {
	      if( initialDecisions.at(i_l1t).second ) l1_doublemu0_eta1p4_os = 1;
	    }
	  if( l1tName == "L1_DoubleMu_10_0_dEta_Max1p8" ) 
	    {
	      if( initialDecisions.at(i_l1t).second ) l1_doublemu_10_0 = 1;
	    }
	  if( l1tName == "L1_DoubleMu_11_4" ) 
	    {
	      if( initialDecisions.at(i_l1t).second ) l1_doublemu_11_4 = 1;
	    }
	}
    }
  else
    {  // data 
      //ESHandle<L1TUtmTriggerMenu> l1GtMenu;
      //iSetup.get<L1TUtmTriggerMenuRcd>().get(l1GtMenu);
      ESHandle<L1TGlobalPrescalesVetos> psAndVetos;
      auto psRcd = iSetup.tryToGet<L1TGlobalPrescalesVetosRcd>();
      if(psRcd) psRcd->get(psAndVetos);
      int columnN= gtUtil_->prescaleColumn();

      for (size_t i_l1t = 0; i_l1t < initialDecisions.size(); i_l1t++) {

	string l1tName = (initialDecisions.at(i_l1t)).first;
	if( l1tName == "NULL") continue;

	if( l1tName == "L1_DoubleMu0" ) {
	  if( initialDecisions.at(i_l1t).second ) l1_doublemu0 = 1;
	  //cout<<(psAndVetos->prescale_table_)[columnN][i_l1t]<<endl;
	}
	if( l1tName == "L1_TripleMu0" ) {
	  if( initialDecisions.at(i_l1t).second ) l1_triplemu0 = 1;
	  prescale_triplemu0 = (psAndVetos->prescale_table_)[columnN][i_l1t];
	}
	if( l1tName == "L1_TripleMu_5_0_0" ) {
	  if( initialDecisions.at(i_l1t).second ) l1_triplemu500 = 1;
	  prescale_triplemu500 = (psAndVetos->prescale_table_)[columnN][i_l1t];
	}
	if( l1tName == "L1_DoubleMu0er1p6_dEta_Max1p8" ) {
	  if( initialDecisions.at(i_l1t).second ) l1_doublemu0_eta1p6 = 1;
	  prescale_doublemu0_eta1p6 = (psAndVetos->prescale_table_)[columnN][i_l1t];
	}
	if( l1tName == "L1_DoubleMu0er1p6_dEta_Max1p8_OS" ) {
	  if( initialDecisions.at(i_l1t).second ) l1_doublemu0_eta1p6_os = 1;
	  prescale_doublemu0_eta1p6_os = (psAndVetos->prescale_table_)[columnN][i_l1t];
	}
	if( l1tName == "L1_DoubleMu0er1p4_dEta_Max1p8_OS" ) {
	  if( initialDecisions.at(i_l1t).second ) l1_doublemu0_eta1p4_os = 1;
	  prescale_doublemu0_eta1p4_os = (psAndVetos->prescale_table_)[columnN][i_l1t];
	}
	if( l1tName == "L1_DoubleMu_10_0_dEta_Max1p8" ) {
	  if( initialDecisions.at(i_l1t).second ) l1_doublemu_10_0 = 1;
	  prescale_doublemu_10_0 = (psAndVetos->prescale_table_)[columnN][i_l1t];
	}
	if( l1tName == "L1_DoubleMu_11_4" ) {
	  if( initialDecisions.at(i_l1t).second ) l1_doublemu_11_4 = 1;
	  prescale_doublemu_11_4 = (psAndVetos->prescale_table_)[columnN][i_l1t];
	}
      }
    } // data
}


void T3MNtuple::fillEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Event_EventNumber = iEvent.id().event();
  Event_RunNumber = iEvent.id().run();
  Event_bunchCrossing = iEvent.bunchCrossing();
  Event_orbitNumber = iEvent.orbitNumber();
  Event_luminosityBlock = iEvent.luminosityBlock();
  Event_isRealData = iEvent.isRealData();

}


// ------------ method called once each job just before starting event loop  ------------
void 
T3MNtuple::beginJob()
{

  std::cout<<" ------------------------>>>>    T3MNtuple begin Job "<<std::endl;
  cnt_=0;

  Service<TFileService> fs;
  h_n3mu = fs->make<TH1F>("n3mu", "", 10, 0, 10);
  h_step = fs->make<TH1F>("step", "", 10, 0, 10);
  tree = fs->make<TTree>("data", "");


  //=============== Event Block ==============
  tree->Branch("Event_EventNumber", &Event_EventNumber);
  tree->Branch("Event_RunNumber", &Event_RunNumber);
  tree->Branch("Event_bunchCrossing", &Event_bunchCrossing);
  tree->Branch("Event_orbitNumber", &Event_orbitNumber);
  tree->Branch("Event_luminosityBlock", &Event_luminosityBlock);
  tree->Branch("Event_isRealData", &Event_isRealData);
  tree->Branch("puN", &puN, "puN/D");
  //tree->Branch("filterbadGlbMuon", &filterbadGlbMuon, "filterbadGlbMuon/D");
  tree->Branch("gen_flavor", &gen_flavor, "gen_flavor/D");
  tree->Branch("nmu_mom", &nmu_mom, "nmu_mom/D");
  tree->Branch("hlt_doublemu4_lmnrt", &hlt_doublemu4_lmnrt, "hlt_doublemu4_lmnrt/D");
  tree->Branch("hlt_doublemu3_tau3mu", &hlt_doublemu3_tau3mu, "hlt_doublemu3_tau3mu/D");
  tree->Branch("l1_triplemu0", &l1_triplemu0, "l1_triplemu0/D");
  tree->Branch("l1_doublemu0", &l1_doublemu0, "l1_doublemu0/D");
  tree->Branch("l1_doublemu0_eta1p6_os", &l1_doublemu0_eta1p6_os, "l1_doublemu0_eta1p6_os/D");
  tree->Branch("l1_doublemu0_eta1p4_os", &l1_doublemu0_eta1p4_os, "l1_doublemu0_eta1p4_os/D");
  tree->Branch("l1_doublemu0_eta1p6", &l1_doublemu0_eta1p6, "l1_doublemu0_eta1p6/D");
  tree->Branch("l1_doublemu_10_0", &l1_doublemu_10_0, "l1_doublemu_10_0/D");
  tree->Branch("l1_doublemu_11_4", &l1_doublemu_11_4, "l1_doublemu_11_4/D");
  tree->Branch("l1_triplemu500", &l1_triplemu500, "l1_triplemu500/D");
  tree->Branch("prescale_triplemu0", &prescale_triplemu0, "prescale_triplemu0/D");
  tree->Branch("prescale_doublemu0_eta1p6", &prescale_doublemu0_eta1p6, "prescale_doublemu0_eta1p6/D");
  tree->Branch("prescale_doublemu_10_0", &prescale_doublemu_10_0, "prescale_doublemu_10_0/D");
  tree->Branch("prescale_triplemu500", &prescale_triplemu500, "prescale_triplemu500/D");
  tree->Branch("prescale_doublemu0_eta1p6_os", &prescale_doublemu0_eta1p6_os, "prescale_doublemu0_eta1p6_os/D");
  tree->Branch("prescale_doublemu0_eta1p4_os", &prescale_doublemu0_eta1p4_os, "prescale_doublemu0_eta1p4_os/D");
  tree->Branch("prescale_doublemu_11_4", &prescale_doublemu_11_4, "prescale_doublemu_11_4/D");
  tree->Branch("n_reco", &n_reco, "n_reco/I");
  tree->Branch("ifar", &ifar, "ifar/I");
  tree->Branch("ipv_gen", &ipv_gen, "ipv_gen/I");
  tree->Branch("ipv1", &ipv1, "ipv1/I");
  tree->Branch("ipv2", &ipv2, "ipv2/I");
  tree->Branch("pdgid_reco", pdgid_reco, "pdgid_reco[3]/D");
  tree->Branch("momid_reco", momid_reco, "momid_reco[3]/D");
  tree->Branch("vxy_reco", vxy_reco, "vxy_reco[3]/D");
  tree->Branch("n_vtx", &n_vtx, "n_vtx/D");
  tree->Branch("trigmat_reco", trigmat_reco, "trigmat_reco[3]/D");
  tree->Branch("trigmat_new", &trigmat_new, "trigmat_new/D");
  tree->Branch("p_reco", p_reco, "p_reco[3]/D");
  tree->Branch("pt_reco", pt_reco, "pt_reco[3]/D");
  tree->Branch("pt_max", &pt_max, "pt_max/D");
  tree->Branch("pt_min", &pt_min, "pt_min/D");
  tree->Branch("eta_reco", eta_reco, "eta_reco[3]/D");
  tree->Branch("eta_min", &eta_min, "eta_min/D");
  tree->Branch("eta_max", &eta_max, "eta_max/D");
  tree->Branch("phi_reco", phi_reco, "phi_reco[3]/D");
  tree->Branch("p_in", p_in, "p_in[3]/D");
  tree->Branch("eta_in", eta_in, "eta_in[3]/D");
  tree->Branch("phi_in", phi_in, "phi_in[3]/D");
  tree->Branch("p_out", p_out, "p_out[3]/D");
  tree->Branch("eta_out", eta_out, "eta_out[3]/D");
  tree->Branch("phi_out", phi_out, "phi_out[3]/D");
  tree->Branch("p_glb", p_glb, "p_glb[3]/D");
  tree->Branch("eta_glb", eta_glb, "eta_glb[3]/D");
  tree->Branch("phi_glb", phi_glb, "phi_glb[3]/D");
  tree->Branch("pt3mu_reco", &pt3mu_reco, "pt3mu_reco/D");
  tree->Branch("pt2mu_12", &pt2mu_12, "pt2mu_12/D");
  tree->Branch("p3mu_reco", &p3mu_reco, "p3mu_reco/D");
  tree->Branch("eta3mu_reco", &eta3mu_reco, "eta3mu_reco/D");
  tree->Branch("m3mu_reco", &m3mu_reco, "m3mu_reco/D");
  tree->Branch("m3mu_simp", &m3mu_simp, "m3mu_simp/D");
  tree->Branch("m3mu_refit", &m3mu_refit, "m3mu_refit/D");
  tree->Branch("pf_reco", pf_reco, "pf_reco[3]/D");
  tree->Branch("rpcmu_reco", rpcmu_reco, "rpcmu_reco[3]/D");
  tree->Branch("rpchits_reco", rpchits_reco, "rpchits_reco[3]/D");
  tree->Branch("comp2d_reco", comp2d_reco, "comp2d_reco[3]/D");
  tree->Branch("tma_reco", tma_reco, "tma_reco[3]/D");
  tree->Branch("tmost_reco", tmost_reco, "tmost_reco[3]/D");
  tree->Branch("tmosat_reco", tmosat_reco, "tmosat_reco[3]/D");
  tree->Branch("tmlst_reco", tmlst_reco, "tmlst_reco[3]/D");
  tree->Branch("tmlsat_reco", tmlsat_reco, "tmlsat_reco[3]/D");
  tree->Branch("tmlsolpt_reco", tmlsolpt_reco, "tmlsolpt_reco[3]/D");
  tree->Branch("tmlsoblpt_reco", tmlsoblpt_reco, "tmlsoblpt_reco[3]/D");
  tree->Branch("calocomp_reco", calocomp_reco, "calocomp_reco[3]/D");
  tree->Branch("segmcomp_reco", segmcomp_reco, "segmcomp_reco[3]/D");
  tree->Branch("segmcomp_0", segmcomp_0, "segmcomp_0[3]/D");
  tree->Branch("segmcomp_1", segmcomp_1, "segmcomp_1[3]/D");
  tree->Branch("segmcomp_2", segmcomp_2, "segmcomp_2[3]/D");
  tree->Branch("segmcomp_3", segmcomp_3, "segmcomp_3[3]/D");
  //tree->Branch("segmcomp_4", segmcomp_4, "segmcomp_4[3]/D");
  //tree->Branch("segmcomp_5", segmcomp_5, "segmcomp_5[3]/D");
  //tree->Branch("segmcomp_6", segmcomp_6, "segmcomp_6[3]/D");
  //tree->Branch("segmcomp_7", segmcomp_7, "segmcomp_7[3]/D");
  tree->Branch("trkhp_reco", trkhp_reco, "trkhp_reco[3]/D");
  tree->Branch("glbnC_reco", glbnC_reco, "glbnC_reco[3]/D");
  tree->Branch("nOVMH_reco", nOVMH_reco, "nOVMH_reco[3]/D");
  tree->Branch("nOVPH_reco", nOVPH_reco, "nOVPH_reco[3]/D");
  tree->Branch("nOVTH_reco", nOVTH_reco, "nOVTH_reco[3]/D");
  tree->Branch("nOLTH_reco", nOLTH_reco, "nOLTH_reco[3]/D");
  tree->Branch("nOLTHin_reco", nOLTHin_reco, "nOLTHin_reco[3]/D");
  tree->Branch("nOLTHout_reco", nOLTHout_reco, "nOLTHout_reco[3]/D");
  //tree->Branch("nSCH_reco", nSCH_reco, "nSCH_reco[3]/D");
  tree->Branch("iTvF_reco", iTvF_reco, "iTvF_reco[3]/D");
  tree->Branch("tLWM_reco", tLWM_reco, "tLWM_reco[3]/D");
  tree->Branch("pLWM_reco", pLWM_reco, "pLWM_reco[3]/D");
  tree->Branch("ggm_reco", ggm_reco, "ggm_reco[3]/D");
  tree->Branch("tgm_reco", tgm_reco, "tgm_reco[3]/D");
  tree->Branch("pterr_reco", pterr_reco, "pterr_reco[3]/D");
  tree->Branch("Iso_sumPt", Iso_sumPt, "Iso_sumPt[3]/D");
  tree->Branch("Iso_nTr", Iso_nTr, "Iso_nTr[3]/D");
  tree->Branch("Iso_emEt", Iso_emEt, "Iso_emEt[3]/D");
  tree->Branch("Iso_hadEt", Iso_hadEt, "Iso_hadEt[3]/D");
  tree->Branch("Iso_eVE", Iso_eVE, "Iso_eVE[3]/D");
  tree->Branch("Iso_hVE", Iso_hVE, "Iso_hVE[3]/D");
  tree->Branch("nOMS_reco", nOMS_reco, "nOMS_reco[3]/D");
  tree->Branch("nOM_reco", nOM_reco, "nOM_reco[3]/D");
  tree->Branch("mSWVH_reco", mSWVH_reco, "mSWVH_reco[3]/D");
  tree->Branch("cschits_sta1", cschits_sta1, "cschits_sta1[3]/D");
  tree->Branch("cscchi2_sta1", cscchi2_sta1, "cscchi2_sta1[3]/D");
  tree->Branch("cschits_sta2", cschits_sta2, "cschits_sta2[3]/D");
  tree->Branch("cscchi2_sta2", cscchi2_sta2, "cscchi2_sta2[3]/D");
  tree->Branch("cscdxdz_sta1", cscdxdz_sta1, "cscdxdz_sta1[3]/D");
  tree->Branch("cscdxdz_sta2", cscdxdz_sta2, "cscdxdz_sta2[3]/D");
  tree->Branch("cscdydz_sta1", cscdydz_sta1, "cscdydz_sta1[3]/D");
  tree->Branch("cscdydz_sta2", cscdydz_sta2, "cscdydz_sta2[3]/D");
  tree->Branch("cscnsegm_sta1", cscnsegm_sta1, "cscnsegm_sta1[3]/D");
  tree->Branch("cscnsegm_sta2", cscnsegm_sta2, "cscnsegm_sta2[3]/D");
  tree->Branch("nOMS_min", &nOMS_min, "nOMS_min/D");
  tree->Branch("tLWM_min", &tLWM_min, "tLWM_min/D");
  tree->Branch("d0_reco", d0_reco, "d0_reco[3]/D");
  tree->Branch("d0sig_reco", d0sig_reco, "d0sig_reco[3]/D");
  tree->Branch("dz_reco", dz_reco, "dz_reco[3]/D");
  tree->Branch("dzpv_reco", dzpv_reco, "dzpv_reco[3]/D");
  tree->Branch("dr12_reco", &dr12_reco, "dr12_reco/D");
  tree->Branch("dr23_reco", &dr23_reco, "dr23_reco/D");
  tree->Branch("dr31_reco", &dr31_reco, "dr31_reco/D");
  tree->Branch("dr_min", &dr_min, "dr_min/D");
  tree->Branch("drtau_max", &drtau_max, "drtau_max/D");
  tree->Branch("charge_reco", charge_reco, "charge_reco[3]/D");
  tree->Branch("drtau_reco", drtau_reco, "drtau_reco[3]/D");
  tree->Branch("dca12_reco", &dca12_reco, "dca12_reco/D");
  tree->Branch("dca23_reco", &dca23_reco, "dca23_reco/D");
  tree->Branch("dca31_reco", &dca31_reco, "dca31_reco/D");
  tree->Branch("dca_max", &dca_max, "dca_max/D");
  tree->Branch("trkrel_reco", trkrel_reco, "trkrel_reco[3]/D");
  tree->Branch("trkrel_max", &trkrel_max, "trkrel_max/D");
  tree->Branch("trkrel_tau", &trkrel_tau, "trkrel_tau/D");
  tree->Branch("trkrel_tau05", &trkrel_tau05, "trkrel_tau05/D");
  tree->Branch("ntrk_reco", ntrk_reco, "ntrk_reco[3]/D");
  tree->Branch("ntrk_sum", &ntrk_sum, "ntrk_sum/D");
  tree->Branch("ntrk_tau", &ntrk_tau, "ntrk_tau/D");
  tree->Branch("mindca_iso", &mindca_iso, "mindca_iso/D");
  tree->Branch("ntrk_tau05", &ntrk_tau05, "ntrk_tau05/D");
  tree->Branch("ntrk_tau_b", &ntrk_tau_b, "ntrk_tau_b/D");
  tree->Branch("mindca_iso05", &mindca_iso05, "mindca_iso05/D");
  tree->Branch("fv_tC", &fv_tC, "fv_tC/D");
  tree->Branch("fv_nC", &fv_nC, "fv_nC/D");
  tree->Branch("fvwo_tC", fvwo_tC, "fvwo_tC[3]/D");
  tree->Branch("fvwo_nC", fvwo_nC, "fvwo_nC[3]/D");
  tree->Branch("fv_Prob", &fv_Prob, "fv_Prob/D");
  tree->Branch("fv_dxy", &fv_dxy, "fv_dxy/D");
  tree->Branch("fv_dxysig", &fv_dxysig, "fv_dxysig/D");
  tree->Branch("fv_ppdl", &fv_ppdl, "fv_ppdl/D");
  tree->Branch("fv_cosdphi", &fv_cosdphi, "fv_cosdphi/D");
  tree->Branch("fv_d3D", &fv_d3D, "fv_d3D/D");
  tree->Branch("fv_d3Dsig", &fv_d3Dsig, "fv_d3Dsig/D");
  tree->Branch("fv_ppdl3D", &fv_ppdl3D, "fv_ppdl3D/D");
  tree->Branch("fv_cosdphi3D", &fv_cosdphi3D, "fv_cosdphi3D/D");
  tree->Branch("iTnC_reco", iTnC_reco, "iTnC1_reco[3]/D");
  tree->Branch("cLP_reco", cLP_reco, "cLP_reco[3]/D");
  tree->Branch("cLM_reco", cLM_reco, "cLM_reco[3]/D");
  tree->Branch("tKink_reco", tKink_reco, "tKink_reco[3]/D");
  tree->Branch("gKink_reco", gKink_reco, "gKink_reco[3]/D");
  tree->Branch("qprod_reco", qprod_reco, "qprod_reco[3]/D");
  tree->Branch("vmuonhitcomb_reco", vmuonhitcomb_reco, "vmuonhitcomb_reco[3]/D");
  tree->Branch("outerchi2_reco", outerchi2_reco, "outerchi2_reco[3]/D");
  tree->Branch("timeatipinouterr_reco", timeatipinouterr_reco, "timeatipinouterr_reco[3]/D");
  tree->Branch("uSta_reco", uSta_reco, "uSta_reco[3]/D");
  tree->Branch("tRC2_reco", tRC2_reco, "tRC2_reco[3]/D");
  tree->Branch("sRC2_reco", sRC2_reco, "sRC2_reco[3]/D");
  tree->Branch("lDist_reco", lDist_reco, "lDist_reco[3]/D");
  tree->Branch("gDEP_reco", gDEP_reco, "gDEP_reco[3]/D");
  tree->Branch("tMat_reco", tMat_reco, "tMat_reco[3]/D");
  tree->Branch("gTP_reco", gTP_reco, "gTP_reco[3]/D");
  tree->Branch("calem_reco", calem_reco, "calem_reco[3]/D");
  tree->Branch("calemS9_reco", calemS9_reco, "calemS9_reco[3]/D");
  tree->Branch("calemS25_reco", calemS25_reco, "calemS25_reco[3]/D");
  tree->Branch("calhad_reco", calhad_reco, "calhad_reco[3]/D");
  tree->Branch("calhadS9_reco", calhadS9_reco, "calhadS9_reco[3]/D");
  tree->Branch("ddz_12", &ddz_12, "ddz_12/D");
  tree->Branch("calomuon_3", &calomuon_3, "calomuon_3/D");
  tree->Branch("n_sv", &n_sv, "n_sv/I");
  tree->Branch("sv_mass", sv_mass, "sv_mass[n_sv]/D");
  tree->Branch("sv_nmu", sv_nmu, "sv_nmu[n_sv]/D");
  tree->Branch("sv_d3D", sv_d3D, "sv_d3D[n_sv]/D");
  tree->Branch("sv_d3Dsig", sv_d3Dsig, "sv_d3Dsig[n_sv]/D");
  tree->Branch("sv_ppdl3D", sv_ppdl3D, "sv_ppdl3D[n_sv]/D");
  tree->Branch("sv_ntrk", sv_ntrk, "sv_ntrk[n_sv]/D");
  tree->Branch("sv_pt", sv_pt, "sv_pt[n_sv]/D");
  tree->Branch("sv_dz", sv_dz, "sv_dz[n_sv]/D");
  tree->Branch("sv_cosdphi3D", sv_cosdphi3D, "sv_cosdphi3D[n_sv]/D");
  tree->Branch("sv_overlap", sv_overlap, "sv_overlap[n_sv]/D");
  tree->Branch("pv_nmu", &pv_nmu, "pv_nmu/D");
  tree->Branch("pv1_tC", &pv1_tC, "pv1_tC/D");
  tree->Branch("pv1_nC", &pv1_nC, "pv1_nC/D");
  tree->Branch("pv2_tC", &pv2_tC, "pv2_tC/D");
  tree->Branch("pv2_nC", &pv2_nC, "pv2_nC/D");
  tree->Branch("ntrk0p1", &ntrk0p1, "ntrk0p1/D");
  tree->Branch("ntrk0p2", &ntrk0p2, "ntrk0p2/D");
  tree->Branch("ntrk0p5", &ntrk0p5, "ntrk0p5/D");
  tree->Branch("maxdxy_pv0", &maxdxy_pv0, "maxdxy_pv0/D");
  tree->Branch("njet20", &njet20, "njet20/I");
  tree->Branch("jet_pt", jet_pt, "jet_pt[njet20]/D");
  tree->Branch("jet_overlap", jet_overlap, "jet_overlap[njet20]/D");
  tree->Branch("btagcvsb", btagcvsb, "btagcvsb[njet20]/D");
  tree->Branch("btagcsv", btagcsv, "btagcsv[njet20]/D");
  tree->Branch("btagmva", btagmva, "btagmva[njet20]/D");
  tree->Branch("m2mu_ss", &m2mu_ss, "m2mu_ss/D");
  tree->Branch("m2mu_os1", &m2mu_os1, "m2mu_os1/D");
  tree->Branch("m2mu_os2", &m2mu_os2, "m2mu_os2/D");
  tree->Branch("m2mu_12", &m2mu_12, "m2mu_12/D");
  tree->Branch("m2mu_23", &m2mu_23, "m2mu_23/D");
  tree->Branch("m2mu_31", &m2mu_31, "m2mu_31/D");
  tree->Branch("m2mu_max", &m2mu_max, "m2mu_max/D");
  tree->Branch("m2mu_min", &m2mu_min, "m2mu_min/D");

  tree->Branch("Track_p4", &Track_p4);
  tree->Branch("Track_normalizedChi2", &Track_normalizedChi2);
  tree->Branch("Track_numberOfValidHits", &Track_numberOfValidHits);
  tree->Branch("Track_charge", &Track_charge);
  tree->Branch("Track_dxy", &Track_dxy);
  tree->Branch("Track_dz", &Track_dz);
  tree->Branch("Track_poca", &Track_poca);
  tree->Branch("Track_dxyError", &Track_dxyError);
  tree->Branch("Track_dzError", &Track_dzError);
  //refitter_.setServices(iSetup);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
T3MNtuple::endJob() 
{
  std::cout << " No Of event processed: " << cnt_ << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  T3MNtuple::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  T3MNtuple::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  T3MNtuple::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  T3MNtuple::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
T3MNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}




void T3MNtuple::ClearEvent() {


  puN = 0;

  l1_doublemu0 = 0; l1_triplemu0 = 0; l1_triplemu500 = 0;
  l1_doublemu_10_0 = 0; l1_doublemu_11_4 = 0;
  l1_doublemu0_eta1p6 = 0; l1_doublemu0_eta1p6_os = 0; l1_doublemu0_eta1p4_os = 0;
  prescale_triplemu0 = 0; prescale_doublemu_10_0 = 0; prescale_doublemu0_eta1p6 = 0;
  prescale_triplemu500 = 0; prescale_doublemu_11_4 = 0; prescale_doublemu0_eta1p6_os = 0;
  prescale_doublemu0_eta1p4_os = 0;

  Track_p4.clear();
  Track_normalizedChi2.clear();
  Track_numberOfValidHits.clear();
  Track_charge.clear();
  Track_dxy.clear();
  Track_dz.clear();
  Track_poca.clear();
  Track_dxyError.clear();
  Track_dzError.clear();

}

//define this as a plug-in
DEFINE_FWK_MODULE(T3MNtuple);
