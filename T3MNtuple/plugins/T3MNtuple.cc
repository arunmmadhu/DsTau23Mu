#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
//Simple Fits

#include "DsTau23Mu/T3MNtuple/interface/TrackParticle.h"
#include "DsTau23Mu/T3MNtuple/interface/ParticleBuilder.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//
double T3MNtuple::MuonPtCut_(-1.);
double T3MNtuple::MuonEtaCut_(999);
double T3MNtuple::ElectronPtCut_(-1.);
double T3MNtuple::ElectronEtaCut_(999);
double T3MNtuple::TrackPtCut_(-1.);
double T3MNtuple::TrackEtaCut_(999);
double T3MNtuple::phimassmin_(1.7);
double T3MNtuple::phimassmax_(3.0);


//
// constructors and destructor
//
T3MNtuple::T3MNtuple(const edm::ParameterSet& iConfig):
   TriggerMuonMatchingdr_(iConfig.getUntrackedParameter("TriggerMuonMatchingdr", (double) 0.3)),
   WhatData_(iConfig.getUntrackedParameter<string>("WhatData","2017")),
   btagCvsBToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("btagsCvsB"))),
   btagCSVToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("btagsCSV"))),
   btagMVAToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("btagsMVA"))),
   //btagDeepCSVToken_(consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("btagDeepCSV"))), ## absent in AOD
   vtxToken_(consumes<VertexCollection>(iConfig.getParameter<InputTag>("pvs"))),
   //ksToken_(consumes<CompositeVertexPtrCollection>(iConfig.getParameter<InputTag>("kshort")))
   //conversionToken_(consumes<CompositeiCandidateCollection>(iConfig.getParameter<InputTag>("conversions")))
   trackToken_(consumes<TrackCollection>(iConfig.getParameter<InputTag>("trks"))),
   triggerToken_(consumes<TriggerResults>(iConfig.getParameter<InputTag>("triggerBitsH"))),
   trigeventToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<InputTag>("triggerSummary"))),
   triggerObjectToken_(consumes<vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<InputTag>("triggerObjects"))),
   algToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<InputTag>("AlgInputTag"))),
   bsToken_(consumes<BeamSpot>(iConfig.getParameter<InputTag>("beamSpotHandle"))),
   puToken_(consumes<vector<PileupSummaryInfo> >(iConfig.getParameter<InputTag>("pileupSummary"))),
   genToken_(consumes<GenParticleCollection>(iConfig.getParameter<InputTag>("genParticles"))),
   patMuonToken_(consumes<vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("pat_muons"))),
   patElectronsToken_(consumes<vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("pat_electrons"))),
   pat_met_puppi_(consumes<vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("met_puppi"))),
   patPhotonToken_(consumes<vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("pat_phos"))),
   compositeSVToken_(consumes<vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<InputTag>("composite_svs"))),
   recoMuonToken_(consumes<vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("reco_muons"))),
   recoPhotonToken_(consumes<vector<reco::Photon>>(iConfig.getParameter<edm::InputTag>("reco_phos"))),
   svToken_(consumes<VertexCollection>(iConfig.getParameter<InputTag>("reco_svs"))),
   TauCandidateToken_(consumes<pat::TauCollection>(edm::InputTag("slimmedTausPlusDeepTau"))),
   pfTauToken_(consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer", "", "PAT"))),
   //pfTauToken_(consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer"))),
   //pfTauToken_(consumes<reco::PFTauCollection>(iConfig.getParameter<InputTag>("PFTauTag"))),
   dmfNewToken_(consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs", "", "PAT"))),
   goodPVToken_(consumes<vector<Vertex>>(iConfig.getParameter<InputTag>("goodPVTag"))),
   thePFCandToken_(consumes<edm::View<pat::PackedCandidate>>        (iConfig.getParameter<edm::InputTag>("thePFCandTag"))),
   tracks_Token_(consumes<edm::View<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"))),
   sampleType_(iConfig.getUntrackedParameter<string>("DataMCType",""))
   
{


   gtUtil_ = new L1TGlobalUtil(iConfig, consumesCollector(), *this, algInputTag_, algInputTag_);
   doMC_ = iConfig.getParameter<bool>("doMC");
   doFullMC_ = iConfig.getParameter<bool>("doFullMC");
   wideSB_ = iConfig.getParameter<bool>("wideSB");
   do2mu_ = iConfig.getParameter<bool>("do2mu");
   passhlt_ = iConfig.getParameter<bool>("passhlt");
   mid_ = iConfig.getParameter<int>("mid");
   doTracks_ = iConfig.getParameter<bool>("doTracks");
   doMuons_ = iConfig.getParameter<bool>("doMuons");
   doTaus_ = iConfig.getParameter<bool>("doTaus");
   doElectrons_ = iConfig.getParameter<bool>("doElectrons");
   do3mutuple_ = iConfig.getParameter<bool>("do3mutuple");
   doL1_ = iConfig.getParameter<bool>("doL1");
   doBJets_ = iConfig.getParameter<bool>("doBJets");
   doPhotons_ = iConfig.getParameter<bool>("doPhotons");
   doThreeMuons_=  iConfig.getParameter<bool>("doThreeMuons");
   doTwoMuonsAndTrack_= iConfig.getParameter<bool>("doTwoMuonsAndTrack");
   MuonPtCut_ = iConfig.getParameter<double>("MuonPtCut"); //default: 1.0
   MuonEtaCut_ = iConfig.getParameter<double>("MuonEtaCut"); //default: 2.4
   ElectronPtCut_ = iConfig.getParameter<double>("ElectronPtCut"); //default: 1.0
   ElectronEtaCut_ = iConfig.getParameter<double>("ElectronEtaCut"); //default: 2.4
   TrackPtCut_ = iConfig.getParameter<double>("TrackPtCut"); //default: 1.0
   TrackEtaCut_ = iConfig.getParameter<double>("TrackEtaCut"); //default: 2.4
   miniAODRun_ = iConfig.getParameter<bool>("miniAODRun"); // true if running over MiniAOD

   DataMCType DMT;
   Event_DataMC_Type=DMT.GetType(sampleType_);    

   DEBUG = false;
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
   if(track.pt()>TrackPtCut_){
      if(fabs(track.eta())<TrackEtaCut_){return true;
         if(track.hitPattern().trackerLayersWithMeasurement()>5){
            if(track.hitPattern().pixelLayersWithMeasurement()>1) return true;
         }
      }
   }
   return false;
}


std::vector<int> T3MNtuple::SortByPt(std::vector<TLorentzVector> invec){
   double pt1=invec.at(0).Pt();
   double pt2=invec.at(1).Pt();
   double pt3=invec.at(2).Pt();

   std::vector<int> out;
   int i1,i2,i3;


   if(pt1>pt2)
   {
      if(pt2>pt3)
      {
         i1=0; i2 = 1; i3 = 2;
      }
      else if(pt1>pt3)
      {
         i1=0; i2 = 2; i3 = 1;
      }
      else
      {
         i1=2; i2 = 0; i3 = 1;
      }
   }
   else
   {
      if(pt1>pt3)
      {
         i1=1; i2 = 0; i3 = 2;
      }
      else if(pt2>pt3)
      {
         i1=1; i2 = 2; i3 = 0;
      }
      else
      {
         i1=2; i2 = 1; i3 = 0;
      }
   }

   out.push_back(i1);  out.push_back(i2);  out.push_back(i3);

   return out;
}



bool T3MNtuple::AcceptedMuon(reco::MuonRef RefMuon) {
   if((RefMuon->pt() > MuonPtCut_) && (abs(RefMuon->eta()) < MuonEtaCut_)){
      //    if(   /*RefMuon->isPFMuon() &&*/  ( RefMuon->isGlobalMuon() || RefMuon->isTrackerMuon()))   return true;
      //    if(   /*RefMuon->isPFMuon() &&*/  ( RefMuon->isGlobalMuon() || RefMuon->isTrackerMuon()) and 
      //	  RefMuon->innerTrack().isNonnull() && RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0)   return true;

      //    if( RefMuon->innerTrack().isNonnull() && RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0)   return true;


      if(RefMuon->innerTrack().isNonnull() && RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && 
            ( RefMuon->isGlobalMuon() || RefMuon->isTrackerMuon()))      return true;


   }
   return false;
}

bool T3MNtuple::AcceptedMuon(pat::MuonRef RefMuon) {
   if((RefMuon->pt() > MuonPtCut_) && (abs(RefMuon->eta()) < MuonEtaCut_)){
      //    if(   /*RefMuon->isPFMuon() &&*/  ( RefMuon->isGlobalMuon() || RefMuon->isTrackerMuon()))   return true;
      //    if(   /*RefMuon->isPFMuon() &&*/  ( RefMuon->isGlobalMuon() || RefMuon->isTrackerMuon()) and 
      //	  RefMuon->innerTrack().isNonnull() && RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0)   return true;

      //    if( RefMuon->innerTrack().isNonnull() && RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0)   return true;


      if(RefMuon->innerTrack().isNonnull() && RefMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && 
            ( RefMuon->isGlobalMuon() || RefMuon->isTrackerMuon()))      return true;


   }
   return false;
}

bool T3MNtuple::getTrackMatch(const Handle<std::vector<reco::Track> > & trackCollection, reco::TrackRef & refTrack, int &match) {
   match = -1;

   for (unsigned int iTrack = 0; iTrack < trackCollection->size(); iTrack++) {
      const reco::TrackRef Track(trackCollection, iTrack);
      if (refTrack == Track) {
         match = iTrack;
         return true;
      }
   }

   return false;
}


bool T3MNtuple::SkipThisParticle(const reco::GenParticle &GenPar){

   int id = abs(GenPar.pdgId());
   if(id == 21 || id == 1 || id == 2 || id ==3 || id == 4 || id == 5 || id == 6) return true;
   if(id == 22 && GenPar.p4().Pt() < 0.1 ) return true;
   return false;


}

bool T3MNtuple::isGoodGenParticle(const reco::GenParticle &GenPar){

   int id = abs(GenPar.pdgId());
   if (id == PDGInfo::Ds_plus) return true;
   if (id == PDGInfo::B_plus) return true;
   if (id == PDGInfo::B_0) return true;
   return false;
}

// ------------ method called for each event  ------------
   void
T3MNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 
   if (DEBUG)   std::cout<<" ========================  new event =============== "<< std::endl;
   cnt_++;
   ClearEvent();

   // Clear handles

   recoMuonCollection.clear();
   patMuonCollection.clear();
   btagsCvsB.clear();
   btagsCSV.clear(); 
   btagsMVA.clear();
   btagDeepCSV.clear();
   pvs.clear();
   svs.clear();
   compositeSVs.clear();
   recoPhotons.clear();
   patPhotons.clear();
   trackCollection.clear();
   triggerBitsH.clear();
   triggerSummary.clear();
   triggerObjects.clear();
   beamSpotHandle.clear();
   puInfo.clear();
   genParticles.clear();
   tauHandle.clear();
   pfTaus.clear();
   vertexs.clear();
   pfCandHandle.clear();
   tracksHandle.clear();

   edm::TriggerNames triggerNames;  

   fillEventInfo(iEvent, iSetup);
   fillMET(iEvent, iSetup);

   if(doMC_){
      if ( iEvent.getByToken(genToken_, genParticles) &&
	   iEvent.getByToken(puToken_, puInfo) ) fillMCTruth(iEvent, iSetup, genParticles, puInfo);
      else {
         if (!iEvent.getByToken(genToken_, genParticles)) edm::LogError("") << "[T3MNtuple]: GEN collection does not exist!";
         if (!iEvent.getByToken(puToken_, puInfo)) edm::LogError("") << "[T3MNtuple]: PileUp Info does not exist!";
      }
   }

   if ( iEvent.getByToken(triggerToken_, triggerBitsH) ) triggerNames = iEvent.triggerNames( *triggerBitsH );
   else edm::LogError("") << "[T3MNtuple]: Trigger Bits collection does not exist!";
   if ( miniAODRun_ && !iEvent.getByToken(triggerObjectToken_, triggerObjects)) edm::LogError("") << "[T3MNtuple]: Running plugin on MiniAOD; Trigger Object collection does not exist!";
   if ( !miniAODRun_ && !iEvent.getByToken(trigeventToken_, triggerSummary)) edm::LogError("") << "[T3MNtuple]: Running plugin on AOD; Trigger Summary does not exist!";

   if(doL1_ && triggerBitsH.isValid()){
     if (miniAODRun_ && triggerObjects.isValid())
       //       std::cout<<"  ---- new event  "<< std::endl;
       fillTrigger(iEvent, iSetup, triggerBitsH, triggerSummary, triggerObjects, triggerNames);
     if (!miniAODRun_ && triggerSummary.isValid())
       fillTrigger(iEvent, iSetup, triggerBitsH, triggerSummary, triggerObjects, triggerNames);
   }





   if(doThreeMuons_) {
      if (DEBUG) cout<<" ------------- making three-muon candidates ------------- "<<endl;
      if ( iEvent.getByToken(trackToken_, trackCollection) && iEvent.getByToken(bsToken_, beamSpotHandle)) {
         if (miniAODRun_){
            if (iEvent.getByToken(patMuonToken_, patMuonCollection) && triggerObjects.isValid()){
               if (DEBUG) cout<<"miniaod 3 muons"<<endl;
               Event_nsignal_candidates = fillThreeMuons(iEvent, iSetup, patMuonCollection, trackCollection, beamSpotHandle, triggerObjects, triggerNames);}
            else  edm::LogError("") << "[T3MNtuple]: Running the plugin on MINIAOD; PAT Muon collection does not exist!";
         }
         else {
            if (iEvent.getByToken(recoMuonToken_, recoMuonCollection) && triggerSummary.isValid() && triggerBitsH.isValid())
               Event_nsignal_candidates = fillThreeMuons(iEvent, iSetup, recoMuonCollection, trackCollection, beamSpotHandle, triggerSummary);
            else  edm::LogError("") << "[T3MNtuple]: Running the plugin on AOD; RECO Muon collection does not exist!";
         }
      }
      else if (!trackCollection.isValid()) edm::LogError("") << "[T3MNtuple]: Track collection does not exist!";   
      else if (!beamSpotHandle.isValid()) edm::LogError("") << "[T3MNtuple]:  BeamSpot does not exist!";   
   }

   if(doTwoMuonsAndTrack_) {
      if (DEBUG) cout<<" ------------- making two-muons-and-track candidates ------------- "<<endl;
      if ( iEvent.getByToken(trackToken_, trackCollection) && iEvent.getByToken(bsToken_, beamSpotHandle))
      {
         if (miniAODRun_){
            if (iEvent.getByToken(patMuonToken_, patMuonCollection) && triggerObjects.isValid() && triggerBitsH.isValid())
               Event_ndsphipi_candidate = fillTwoMuonsAndTracks(iEvent, iSetup, patMuonCollection, trackCollection, beamSpotHandle, triggerObjects, triggerNames);
            else  edm::LogError("") << "[T3MNtuple]: PAT Muon collection does not exist!";
         }
         else {
            if (iEvent.getByToken(recoMuonToken_, recoMuonCollection) && triggerSummary.isValid())
               Event_ndsphipi_candidate = fillTwoMuonsAndTracks(iEvent, iSetup, recoMuonCollection, trackCollection, beamSpotHandle, triggerSummary);
            else  edm::LogError("") << "[T3MNtuple]: RECO Muon collection does not exist!";
         }
      }
      else if (!trackCollection.isValid()) edm::LogError("") << "[T3MNtuple]: Track collection does not exist!";       
      else if (!beamSpotHandle.isValid()) edm::LogError("") << "[T3MNtuple]:  BeamSpot does not exist!";   
   }

   if (DEBUG) cout<<"Number of signal candidates = "<<Event_nsignal_candidates<<endl;
   if (DEBUG) cout<<"Number of dsphipi candidates = "<<Event_ndsphipi_candidate<<endl;

   // Fill the output tree only if a candidate is found 
   if(Event_nsignal_candidates!=0 or Event_ndsphipi_candidate!=0){
   //if(Event_nsignal_candidates!=0 or Event_ndsphipi_candidate!=0 or doTaus_){

      MC_isReco=1;
      
      if (!iEvent.getByToken(vtxToken_, pvs)) edm::LogError("") << "[T3MNtuple]: Primary Vertex collection does not exist!";
      if (!miniAODRun_ && !iEvent.getByToken(svToken_, svs)) edm::LogError("") << "[T3MNtuple]: Secondary Vertex collection does not exist!";
      if (miniAODRun_ && !iEvent.getByToken(compositeSVToken_, compositeSVs)) edm::LogError("") << "[T3MNtuple]: Composite Secondary Vertex collection does not exist!";

      if  (beamSpotHandle.isValid() && pvs.isValid() ){
         if (miniAODRun_ && patMuonCollection.isValid() && compositeSVs.isValid()){
            if (DEBUG) cout<<"Filling Vertex information using PAT muons.";
            fillVertices(iEvent, iSetup, patMuonCollection, trackCollection, pvs, compositeSVs, beamSpotHandle);
         }
         else if (!miniAODRun_ && recoMuonCollection.isValid() && svs.isValid()){
            if (DEBUG) cout<<"Filling Vertex information using RECO muons.";
            fillVertices(iEvent, iSetup, recoMuonCollection, trackCollection, pvs, svs, beamSpotHandle);
         }
      }

      if(doBJets_ && !miniAODRun_){
         if ( !miniAODRun_ && iEvent.getByToken(btagCvsBToken_, btagsCvsB) &&
               iEvent.getByToken(btagCSVToken_, btagsCSV) &&
               iEvent.getByToken(btagMVAToken_, btagsMVA) ) {
            if (DEBUG) cout<<"Filling BTagJets."<<endl;
            fillBTagJets(iEvent, iSetup, btagsCvsB, btagsCSV, btagsMVA);
         }
         else{
            if (!btagsCvsB.isValid()) edm::LogError("") << "[T3MNtuple]: BTagCvsB collection does not exist!";
            if (!btagsCSV.isValid()) edm::LogError("") << "[T3MNtuple]: BTagCSV collection does not exist!";
            if (!btagsMVA.isValid()) edm::LogError("") << "[T3MNtuple]: BTagMVA collection does not exist!";
         }
      }

      if(doPhotons_){  
         if (miniAODRun_ && iEvent.getByToken(patPhotonToken_, patPhotons)) fillPhotons(iEvent, iSetup, patPhotons);
         else if (!miniAODRun_ && iEvent.getByToken(recoPhotonToken_, recoPhotons)) fillPhotons(iEvent, iSetup, recoPhotons);
         else edm::LogError("") << "[T3Mntuple]: Photon collection does not exist!";
      }

      if(doMuons_)
	{      
         if (trackCollection.isValid() && pvs.isValid())
	   {
            if (miniAODRun_ && patMuonCollection.isValid())
	      {
		if (DEBUG) cout<<"Filling PAT muons."<<endl;
		fillMuons(iEvent, iSetup, patMuonCollection, trackCollection, pvs);
	      }
            else if (!miniAODRun_ && recoMuonCollection.isValid()){
	      if (DEBUG) cout<<"Filling RECO muons."<<endl;
	      fillMuons(iEvent, iSetup, recoMuonCollection, trackCollection, pvs);
            }
	   }
	}
      if(doElectrons_)
	{
          iEvent.getByToken(goodPVToken_ , vertexs);
          iEvent.getByToken(patElectronsToken_, patElectronCollection);
          iEvent.getByToken(tracks_Token_, tracksHandle);
	  if (tracksHandle.isValid() && vertexs.isValid())
	    {	
	      if (miniAODRun_)  // && tauHandle.isValid()
		{
		  fillElectrons(iEvent, iSetup, trackCollection, pvs, beamSpotHandle, patElectronCollection, vertexs);
		}
	    }
	}
   
     
      if(doTaus_)
	{      
	  //	  cout<<" Getting tauHandle: "<< 
	  //	    iEvent.getByToken(TauCandidateToken_, tauHandle) <<" Getting pfTaus: "<< 
	  //	    iEvent.getByToken(pfTauToken_, pfTaus) <<" Getting dmfNew: "<< 
	  //	    iEvent.getByToken(dmfNewToken_, dmfNew) <<" Getting vertexs: "<< 
	  //	    iEvent.getByToken(goodPVToken_ , vertexs) <<" Getting pfCandHandle: "<< 
	  //	    iEvent.getByToken(thePFCandToken_, pfCandHandle) <<" Getting tracksHandle: "<< iEvent.getByToken(tracks_Token_, tracksHandle) <<endl;

	  iEvent.getByToken(TauCandidateToken_, tauHandle);
	  iEvent.getByToken(goodPVToken_ , vertexs);
	  iEvent.getByToken(thePFCandToken_, pfCandHandle);
	  iEvent.getByToken(tracks_Token_, tracksHandle);
	  if (tracksHandle.isValid() && vertexs.isValid())
	    {	
	      if (miniAODRun_)  // && tauHandle.isValid()
		{
		  fillTaus(iEvent, iSetup, trackCollection, pvs, beamSpotHandle, tauHandle, vertexs, pfCandHandle, tracksHandle);
		}
	    }
	}
      
      
         

      if(doTracks_){
         if (trackCollection.isValid()) fillTracks(iEvent, iSetup, trackCollection);
         else edm::LogError("") << "Track collection does not exist!";
      }

   }
   output_tree->Fill();

}


std::vector<const reco::GenParticle* > T3MNtuple::TauDecayProducts(const reco::GenParticle *Tau){
   std::vector<const reco::GenParticle* > out;
   unsigned int pdgid=abs(Tau->pdgId());
   if(pdgid==PDGInfo::tau_minus){ // check that it is a tau
      out.push_back(Tau);
      for (unsigned int i=0; i< Tau->numberOfDaughters(); i++){
         const reco::Candidate *dau=Tau->daughter(i);
         out.push_back(static_cast<const reco::GenParticle*>(dau));
      }
   }
   return out;
}

void T3MNtuple::matchLegCombinations(vector<TLorentzVector>& particleP4, vector<TLorentzVector>& trigObjectP4, double drmax, vector<float>& match)
{

   int SIZE = trigObjectP4.size();

   if (SIZE<3) return;
   if (DEBUG) cout<<"Particle size = "<<particleP4.size()<<endl;

   double dpt[3][SIZE];
   double dR[3][SIZE];
   
   for (int i=0; i<3; i++){
      for (int j=0; j<SIZE; j++){
         double p = fabs(particleP4[i].Pt()-trigObjectP4[j].Pt())/particleP4[i].Pt();
         double r = particleP4[i].DeltaR(trigObjectP4[j]);
         dpt[i][j] = p;
         dR[i][j] = r;
      }
   }

   // Find if any combination matches the trigger matching criteria
   for (int i=0; i<SIZE; i++){
      for (int j=0; j<SIZE; j++){
         if (i==j) continue;
         for (int k=0; k<SIZE; k++){
            if (k==i || k==j) continue;
            if ( dpt[0][i]<0.1 && dpt[1][j]<0.1 && dpt[2][k]<0.1 &&
                 dR[0][i]<drmax && dR[1][j]<drmax && dR[2][k]<drmax ) {
               match.at(0) = dR[0][i];
               match.at(1) = dR[1][j];
               match.at(2) = dR[2][k];
            }
         }
      }
   }

   if (DEBUG) cout<<"Trigger matching complete!"<<endl; 
}


void T3MNtuple::TriggerMatch(const Handle<trigger::TriggerEvent> &triggerSummary, 
      vector<TLorentzVector>& particleP4, double drmax, vector<float> &match) {
   
   std::vector<trigger::TriggerObject> trigobjs = triggerSummary->getObjects();
   edm::InputTag MuonFilterTag;

   if(WhatData_=="2017")  MuonFilterTag = edm::InputTag("hltTau3muTkVertexFilter", "", "HLT"); 
   if(WhatData_=="2018")  MuonFilterTag = edm::InputTag("hltdstau3muDisplaced3muFltr", "", "HLT"); 
   //  if(WhatData_=="Parked")  edm::InputTag MuonFilterTag = edm::InputTag("hltdstau3muDisplaced3muFltr", "", "HLT");
   
   vector<TLorentzVector> trigObjectP4;
   
   size_t MuonFilterIndex = (*triggerSummary).filterIndex(MuonFilterTag); 
   if(MuonFilterIndex < (*triggerSummary).sizeFilters()) {
      const trigger::Keys &KEYS = (*triggerSummary).filterKeys(MuonFilterIndex);
      for (unsigned int ipart = 0; ipart < KEYS.size(); ipart++) {

         double tmp_p = trigobjs.at(KEYS.at(ipart)).p();

         TLorentzVector tmpP4(trigobjs.at(KEYS.at(ipart)).px(),
                              trigobjs.at(KEYS.at(ipart)).py(),
                              trigobjs.at(KEYS.at(ipart)).pz(),
                              sqrt(pow(tmp_p,2) + pow(PDGInfo::mu_mass(),2)) );

         trigObjectP4.push_back(tmpP4);
      }
   }
   
   if (DEBUG) cout<<"Matching Trigger legs to tracks"<<endl; 
   matchLegCombinations(particleP4, trigObjectP4, drmax, match);
}

void T3MNtuple::TriggerMatch(const edm::Event& iEvent, const Handle<vector<pat::TriggerObjectStandAlone>>& triggerObjects, const TriggerNames& triggerNames, 
      vector<TLorentzVector>& particleP4, double drmax, vector<float> &match){

   std::string MuonFilterTag;

   if(WhatData_=="2017")  MuonFilterTag = "hltTau3muTkVertexFilter"; 
   if(WhatData_=="2018")  MuonFilterTag = "hltdstau3muDisplaced3muFltr"; 

   vector<TLorentzVector> trigObjectP4;
   
   for (pat::TriggerObjectStandAlone to: *triggerObjects){

      if (DEBUG) cout<<"unpacking names"<<endl;
      to.unpackPathNames(triggerNames);
      trigger::size_type nFilters = to.filterLabels().size();
      for (trigger::size_type iFilter=0; iFilter<nFilters; iFilter++){
         std::string filterName = to.filterLabels()[iFilter];

         if (filterName.compare(MuonFilterTag)!=0) continue;

         double tmp_p = to.p();
         TLorentzVector tmpP4(to.px(), to.py(), to.pz(), sqrt(pow(tmp_p,2) + pow(PDGInfo::mu_mass(),2)) );
         trigObjectP4.push_back(tmpP4);
      }
   }

   if (DEBUG) cout<<"Matching Trigger legs to tracks"<<endl; 
   matchLegCombinations(particleP4, trigObjectP4, drmax, match);
}

void T3MNtuple::fillEventInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   Event_EventNumber = iEvent.id().event();
   Event_RunNumber = iEvent.id().run();
   Event_bunchCrossing = iEvent.bunchCrossing();
   Event_orbitNumber = iEvent.orbitNumber();
   Event_luminosityBlock = iEvent.luminosityBlock();
   Event_isRealData = iEvent.isRealData();

   DataMCType DMT;
   Event_DataMC_Type = DMT.GetType();
   if (Event_isRealData) {
      Event_DataMC_Type = DataMCType::Data;
   }
} 



void T3MNtuple::fillMET(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  edm::Handle<std::vector<pat::MET> > METpuppi;
  iEvent.getByToken(pat_met_puppi_, METpuppi);

  Event_METEt = METpuppi->front().et();
  Event_METPhi = METpuppi->front().phi();
  Event_METXX = METpuppi->front().getSignificanceMatrix()(0,0);
  Event_METXY = METpuppi->front().getSignificanceMatrix()(0,1);
  Event_METYY = METpuppi->front().getSignificanceMatrix()(1,1);

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


   output_tree = fs->make<TTree>("t3mtree", "");

   //=============== Event Block ==============
   output_tree->Branch("Event_EventNumber", &Event_EventNumber);
   output_tree->Branch("Event_RunNumber", &Event_RunNumber);
   output_tree->Branch("Event_bunchCrossing", &Event_bunchCrossing);
   output_tree->Branch("Event_orbitNumber", &Event_orbitNumber);
   output_tree->Branch("Event_luminosityBlock", &Event_luminosityBlock);
   output_tree->Branch("Event_isRealData", &Event_isRealData);
   output_tree->Branch("Event_nsignal_candidates", &Event_nsignal_candidates);
   output_tree->Branch("Event_ndsphipi_candidate", &Event_ndsphipi_candidate);
   output_tree->Branch("Event_DataMC_Type" ,&Event_DataMC_Type);
   output_tree->Branch("Event_METEt" ,&Event_METEt);
   output_tree->Branch("Event_METPhi" ,&Event_METPhi);
   output_tree->Branch("Event_METXX" ,&Event_METXX);
   output_tree->Branch("Event_METXY" ,&Event_METXY);
   output_tree->Branch("Event_METYY" ,&Event_METYY);


   output_tree->Branch("puN", &puN, "puN/D");

   output_tree->Branch("Track_p4", &Track_p4);
   output_tree->Branch("Track_normalizedChi2", &Track_normalizedChi2);
   output_tree->Branch("Track_numberOfValidHits", &Track_numberOfValidHits);
   output_tree->Branch("Track_charge", &Track_charge);
   output_tree->Branch("Track_dxy", &Track_dxy);
   output_tree->Branch("Track_dz", &Track_dz);
   output_tree->Branch("Track_poca", &Track_poca);
   output_tree->Branch("Track_dxyError", &Track_dxyError);
   output_tree->Branch("Track_dzError", &Track_dzError);


   //=============  Electrons Block ====
   output_tree->Branch("Electron_p4", &Electron_p4);
   output_tree->Branch("Electron_Charge", &Electron_Charge);
   output_tree->Branch("Electron_puppiNeutralHadronIso", &Electron_puppiNeutralHadronIso);
   output_tree->Branch("Electron_puppiPhotonIso", &Electron_puppiPhotonIso);
   output_tree->Branch("Electron_puppiChargedHadronIso", &Electron_puppiChargedHadronIso);
   output_tree->Branch("Electron_trackIso", &Electron_trackIso);
   output_tree->Branch("Electron_isPF", &Electron_isPF);
   output_tree->Branch("Electron_cutBasedElectronID_Fall17_94X_V2_veto", &Electron_cutBasedElectronID_Fall17_94X_V2_veto);
   output_tree->Branch("Electron_cutBasedElectronID_Fall17_94X_V2_loose", &Electron_cutBasedElectronID_Fall17_94X_V2_loose);
   output_tree->Branch("Electron_cutBasedElectronID_Fall17_94X_V2_medium", &Electron_cutBasedElectronID_Fall17_94X_V2_medium);
   output_tree->Branch("Electron_cutBasedElectronID_Fall17_94X_V2_tight", &Electron_cutBasedElectronID_Fall17_94X_V2_tight);



   //=============  Muon Block ====
   output_tree->Branch("Muon_p4", &Muon_p4);
   output_tree->Branch("Muon_Poca", &Muon_Poca);
   output_tree->Branch("Muon_isGlobalMuon", &Muon_isGlobalMuon);
   output_tree->Branch("Muon_isStandAloneMuon", &Muon_isStandAloneMuon);
   output_tree->Branch("Muon_isTrackerMuon", &Muon_isTrackerMuon);
   output_tree->Branch("Muon_isCaloMuon", &Muon_isCaloMuon);
   output_tree->Branch("Muon_isIsolationValid", &Muon_isIsolationValid);
   output_tree->Branch("Muon_isQualityValid", &Muon_isQualityValid);
   output_tree->Branch("Muon_isTimeValid", &Muon_isTimeValid);
   output_tree->Branch("Muon_expectedNnumberOfMatchedStations", &Muon_expectedNnumberOfMatchedStations);
   output_tree->Branch("Muon_emEt03", &Muon_emEt03);
   output_tree->Branch("Muon_emVetoEt03", &Muon_emVetoEt03);
   output_tree->Branch("Muon_hadEt03", &Muon_hadEt03);
   output_tree->Branch("Muon_hadVetoEt03", &Muon_hadVetoEt03);
   output_tree->Branch("Muon_nJets03", &Muon_nJets03);
   output_tree->Branch("Muon_nTracks03", &Muon_nTracks03);
   output_tree->Branch("Muon_sumPt03", &Muon_sumPt03);
   output_tree->Branch("Muon_trackerVetoPt03", &Muon_trackerVetoPt03);
   output_tree->Branch("Muon_emEt05", &Muon_emEt05);
   output_tree->Branch("Muon_emVetoEt05", &Muon_emVetoEt05);
   output_tree->Branch("Muon_hadEt05", &Muon_hadEt05);
   output_tree->Branch("Muon_hadVetoEt05", &Muon_hadVetoEt05);
   output_tree->Branch("Muon_nJets05", &Muon_nJets05);
   output_tree->Branch("Muon_nTracks05", &Muon_nTracks05);
   output_tree->Branch("Muon_sumPt05", &Muon_sumPt05);
   output_tree->Branch("Muon_timeAtIpInOut", &Muon_timeAtIpInOut);
   output_tree->Branch("Muon_timeAtIpInOutErr", &Muon_timeAtIpInOutErr);


   output_tree->Branch("Muon_trackerVetoPt05", &Muon_trackerVetoPt05);
   output_tree->Branch("Muon_sumChargedHadronPt03", &Muon_sumChargedHadronPt03);
   output_tree->Branch("Muon_sumChargedParticlePt03", &Muon_sumChargedParticlePt03);
   output_tree->Branch("Muon_sumNeutralHadronEt03", &Muon_sumNeutralHadronEt03);
   output_tree->Branch("Muon_sumNeutralHadronEtHighThreshold03", &Muon_sumNeutralHadronEtHighThreshold03);
   output_tree->Branch("Muon_sumPhotonEt03", &Muon_sumPhotonEt03);
   output_tree->Branch("Muon_sumPhotonEtHighThreshold03", &Muon_sumPhotonEtHighThreshold03);
   output_tree->Branch("Muon_sumPUPt03", &Muon_sumPUPt03);
   output_tree->Branch("Muon_sumChargedHadronPt04", &Muon_sumChargedHadronPt04);
   output_tree->Branch("Muon_sumChargedParticlePt04", &Muon_sumChargedParticlePt04);
   output_tree->Branch("Muon_sumNeutralHadronEt04", &Muon_sumNeutralHadronEt04);
   output_tree->Branch("Muon_sumNeutralHadronEtHighThreshold04", &Muon_sumNeutralHadronEtHighThreshold04);
   output_tree->Branch("Muon_sumPhotonEt04", &Muon_sumPhotonEt04);
   output_tree->Branch("Muon_sumPhotonEtHighThreshold04", &Muon_sumPhotonEtHighThreshold04);
   output_tree->Branch("Muon_sumPUPt04", &Muon_sumPUPt04);
   output_tree->Branch("Muon_Track_idx", &Muon_Track_idx);
   output_tree->Branch("Muon_hitPattern_pixelLayerwithMeas", &Muon_hitPattern_pixelLayerwithMeas);
   output_tree->Branch("Muon_numberOfMatchedStations", &Muon_numberOfMatchedStations);
   output_tree->Branch("Muon_normChi2", &Muon_normChi2);
   output_tree->Branch("Muon_hitPattern_numberOfValidMuonHits", &Muon_hitPattern_numberOfValidMuonHits);
   output_tree->Branch("Muon_innerTrack_numberofValidHits", &Muon_innerTrack_numberofValidHits);
   output_tree->Branch("Muon_numberOfMatches", &Muon_numberOfMatches);
   output_tree->Branch("Muon_numberOfChambers", &Muon_numberOfChambers);
   output_tree->Branch("Muon_isPFMuon", &Muon_isPFMuon);
   output_tree->Branch("Muon_isRPCMuon", &Muon_isRPCMuon);
   output_tree->Branch("Muon_numberofValidPixelHits", &Muon_numberofValidPixelHits);
   output_tree->Branch("Muon_trackerLayersWithMeasurement", &Muon_trackerLayersWithMeasurement);

   output_tree->Branch("Muon_combinedQuality_updatedSta",&Muon_combinedQuality_updatedSta);
   output_tree->Branch("Muon_combinedQuality_trkKink",&Muon_combinedQuality_trkKink);
   output_tree->Branch("Muon_combinedQuality_glbKink",&Muon_combinedQuality_glbKink);
   output_tree->Branch("Muon_combinedQuality_trkRelChi2",&Muon_combinedQuality_trkRelChi2);
   output_tree->Branch("Muon_combinedQuality_staRelChi2",&Muon_combinedQuality_staRelChi2);
   output_tree->Branch("Muon_combinedQuality_chi2LocalPosition",&Muon_combinedQuality_chi2LocalPosition);
   output_tree->Branch("Muon_combinedQuality_chi2LocalMomentum",&Muon_combinedQuality_chi2LocalMomentum);
   output_tree->Branch("Muon_combinedQuality_localDistance",&Muon_combinedQuality_localDistance);
   output_tree->Branch("Muon_combinedQuality_globalDeltaEtaPhi",&Muon_combinedQuality_globalDeltaEtaPhi);
   output_tree->Branch("Muon_combinedQuality_tightMatch",&Muon_combinedQuality_tightMatch);
   output_tree->Branch("Muon_combinedQuality_glbTrackProbability",&Muon_combinedQuality_glbTrackProbability);

   output_tree->Branch("Muon_prod_inner_outer_charge",&Muon_prod_inner_outer_charge);
   output_tree->Branch("Muon_outerTrack_p4",&Muon_outerTrack_p4);
   output_tree->Branch("Muon_innerTrack_p4",&Muon_innerTrack_p4);
   output_tree->Branch("Muon_innerTrack_quality",&Muon_innerTrack_quality);
   output_tree->Branch("Muon_ptErrOverPt",&Muon_ptErrOverPt);
   output_tree->Branch("Muon_calEnergy_hadS9",&Muon_calEnergy_hadS9);
   output_tree->Branch("Muon_calEnergy_had",&Muon_calEnergy_had);
   output_tree->Branch("Muon_calEnergy_emS25",&Muon_calEnergy_emS25);
   output_tree->Branch("Muon_calEnergy_emS9",&Muon_calEnergy_emS9);
   output_tree->Branch("Muon_calEnergy_em",&Muon_calEnergy_em);

   output_tree->Branch("Muon_TrackX",&Muon_TrackX);
   output_tree->Branch("Muon_TrackY",&Muon_TrackY);
   output_tree->Branch("Muon_dDxDz",&Muon_dDxDz);
   output_tree->Branch("Muon_dDyDz",&Muon_dDyDz);
   output_tree->Branch("Muon_dX",&Muon_dX);
   output_tree->Branch("Muon_dY",&Muon_dY);
   output_tree->Branch("Muon_pullX",&Muon_pullX);
   output_tree->Branch("Muon_pullY",&Muon_pullY);
   output_tree->Branch("Muon_pullDxDz",&Muon_pullDxDz);
   output_tree->Branch("Muon_pullDyDz",&Muon_pullDyDz);
   output_tree->Branch("numberOfSegments",&numberOfSegments);

   output_tree->Branch("Muon_ptError",&Muon_ptError);
   output_tree->Branch("Muon_phiError",&Muon_phiError);
   output_tree->Branch("Muon_etaError",&Muon_etaError);

   output_tree->Branch("Muon_segmentCompatibility",&Muon_segmentCompatibility);
   output_tree->Branch("Muon_caloCompatibility",&Muon_caloCompatibility);
   output_tree->Branch("Muon_isGoodMuon_TM2DCompatibility",&Muon_isGoodMuon_TM2DCompatibility);

   output_tree->Branch("Muon_innerTrack_validFraction",&Muon_innerTrack_validFraction );
   output_tree->Branch("Muon_innerTrack_pixelLayersWithMeasurement",&Muon_innerTrack_pixelLayersWithMeasurement );
   output_tree->Branch("Muon_innerTrack_numberOfValidTrackerHits",&Muon_innerTrack_numberOfValidTrackerHits );
   output_tree->Branch("Muon_innerTrack_numberOfLostTrackerHits",&Muon_innerTrack_numberOfLostTrackerHits );
   output_tree->Branch("Muon_innerTrack_numberOfLostTrackerInnerHits",&Muon_innerTrack_numberOfLostTrackerInnerHits );
   output_tree->Branch("Muon_innerTrack_numberOfLostTrackerOuterHits",&Muon_innerTrack_numberOfLostTrackerOuterHits );
   output_tree->Branch("Muon_innerTrack_normalizedChi2",&Muon_innerTrack_normalizedChi2 );

   output_tree->Branch("Muon_outerTrack_normalizedChi2",&Muon_outerTrack_normalizedChi2 );
   output_tree->Branch("Muon_outerTrack_muonStationsWithValidHits",&Muon_outerTrack_muonStationsWithValidHits );
   output_tree->Branch("Muon_isGoodMuon_TrackerMuonArbitrated",&Muon_isGoodMuon_TrackerMuonArbitrated );
   output_tree->Branch("Muon_isGoodMuon_TMOneStationTight",&Muon_isGoodMuon_TMOneStationTight );
   output_tree->Branch("Muon_isGoodMuon_TMOneStationAngTight",&Muon_isGoodMuon_TMOneStationAngTight );
   output_tree->Branch("Muon_isGoodMuon_TMLastStationTight",&Muon_isGoodMuon_TMLastStationTight );
   output_tree->Branch("Muon_isGoodMuon_TMLastStationAngTight",&Muon_isGoodMuon_TMLastStationAngTight );
   output_tree->Branch("Muon_isGoodMuon_TMLastStationOptimizedLowPtTight",&Muon_isGoodMuon_TMLastStationOptimizedLowPtTight );
   output_tree->Branch("Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight",&Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight );


   output_tree->Branch("Muon_vmuonhitcomb_reco",&Muon_vmuonhitcomb_reco);
   output_tree->Branch("Muon_rpchits_reco",&Muon_rpchits_reco);

   output_tree->Branch("Muon_ID",&Muon_ID);
   output_tree->Branch("Muon_StandardSelection",&Muon_StandardSelection);

   output_tree->Branch("Muon_charge", &Muon_charge);
   output_tree->Branch("Muon_trackCharge", &Muon_trackCharge);
   output_tree->Branch("Muon_pdgid", &Muon_pdgid);
   output_tree->Branch("Muon_B", &Muon_B);
   output_tree->Branch("Muon_M", &Muon_M);
   output_tree->Branch("Muon_par", &Muon_par);
   output_tree->Branch("Muon_cov", &Muon_cov);





   if(doTaus_){
     output_tree->Branch("Tau_p4", &Tau_p4);
     output_tree->Branch("Tau_charge", &Tau_charge);
     output_tree->Branch("Tau_DecayMode", &Tau_DecayMode);
     output_tree->Branch("Tau_DecayModeFinding", &Tau_DecayModeFinding);
     output_tree->Branch("Tau_NewDecayModeFinding", &Tau_NewDecayModeFinding);

     output_tree->Branch("Tau_byLooseDeepTau2017v2p1VSe", &Tau_byLooseDeepTau2017v2p1VSe);
     output_tree->Branch("Tau_byMediumDeepTau2017v2p1VSe", &Tau_byMediumDeepTau2017v2p1VSe);
     output_tree->Branch("Tau_byTightDeepTau2017v2p1VSe", &Tau_byTightDeepTau2017v2p1VSe);

     output_tree->Branch("Tau_byLooseDeepTau2017v2p1VSmu", &Tau_byLooseDeepTau2017v2p1VSmu);
     output_tree->Branch("Tau_byMediumDeepTau2017v2p1VSmu", &Tau_byMediumDeepTau2017v2p1VSmu);
     output_tree->Branch("Tau_byTightDeepTau2017v2p1VSmu", &Tau_byTightDeepTau2017v2p1VSmu);

     output_tree->Branch("Tau_byLooseDeepTau2017v2p1VSjet", &Tau_byLooseDeepTau2017v2p1VSjet);
     output_tree->Branch("Tau_byMediumDeepTau2017v2p1VSjet", &Tau_byMediumDeepTau2017v2p1VSjet);
     output_tree->Branch("Tau_byTightDeepTau2017v2p1VSjet", &Tau_byTightDeepTau2017v2p1VSjet);


     output_tree->Branch("Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
     output_tree->Branch("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
     output_tree->Branch("Tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &Tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
     output_tree->Branch("Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits);

     output_tree->Branch("Tau_FloatDiscriminants", &Tau_FloatDiscriminants);
     output_tree->Branch("Tau_IntDiscriminants", &Tau_IntDiscriminants);



     output_tree->Branch("Tau_PFTauTrack_p4", &Tau_PFTauTrack_p4);
     output_tree->Branch("Tau_Track_par", &Tau_Track_par);
     output_tree->Branch("Tau_Track_cov", &Tau_Track_cov);

     output_tree->Branch("Tau_Track_Charge", &Tau_Track_Charge);
     output_tree->Branch("Tau_Track_pdgid", &Tau_Track_pdgid);
     output_tree->Branch("Tau_Track_B", &Tau_Track_B);
     output_tree->Branch("Tau_Track_M", &Tau_Track_M);

     output_tree->Branch("Tau_SVPos", &Tau_SVPos);
     output_tree->Branch("Tau_SVCov", &Tau_SVCov);
     output_tree->Branch("Tau_a1_charge", &Tau_a1_charge);
     output_tree->Branch("Tau_a1_pdgid", &Tau_a1_pdgid);
     output_tree->Branch("Tau_a1_B", &Tau_a1_B);
     output_tree->Branch("Tau_a1_M", &Tau_a1_M);
     output_tree->Branch("Tau_a1_lvp", &Tau_a1_lvp);
     output_tree->Branch("Tau_a1_cov", &Tau_a1_cov);
  
   }



   if(doPhotons_){
      output_tree->Branch("Gamma_P4",&Gamma_P4);
      output_tree->Branch("Gamma_hasPixelSeed",&Gamma_hasPixelSeed);
      output_tree->Branch("Gamma_hasConversionTracks",&Gamma_hasConversionTracks);
      output_tree->Branch("Gamma_e1x5",&Gamma_e1x5);
      output_tree->Branch("Gamma_e2x5",&Gamma_e2x5);
      output_tree->Branch("Gamma_e3x3",&Gamma_e3x3);
      output_tree->Branch("Gamma_e5x5",&Gamma_e5x5);
      output_tree->Branch("Gamma_isPFPhoton",&Gamma_isPFPhoton);
   }
   if(doMC_){
      if(doFullMC_){
         output_tree->Branch("MC_p4", &MC_p4);
         output_tree->Branch("MC_vertex", &MC_vertex);
         output_tree->Branch("MC_pdgid", &MC_pdgid);
         output_tree->Branch("MC_charge", &MC_charge);
         output_tree->Branch("MC_midx", &MC_midx);
         output_tree->Branch("MC_childpdgid", &MC_childpdgid);
         output_tree->Branch("MC_childidx", &MC_childidx);
         output_tree->Branch("MC_status", &MC_status);
      }
      output_tree->Branch("MC_isReco", &MC_isReco);
      output_tree->Branch("MCSignalParticle_p4", &MCSignalParticle_p4);
      output_tree->Branch("MCSignalParticle_Vertex", &MCSignalParticle_Vertex);
      output_tree->Branch("MCSignalParticle_pdgid", &MCSignalParticle_pdgid);
      output_tree->Branch("MCSignalParticle_childpdgid", &MCSignalParticle_childpdgid);
      output_tree->Branch("MCSignalParticle_childp4", &MCSignalParticle_childp4);
      output_tree->Branch("MCSignalParticle_Sourcepdgid", &MCSignalParticle_Sourcepdgid);
      output_tree->Branch("MCSignalParticle_Sourcep4", &MCSignalParticle_Sourcep4);
      output_tree->Branch("MCSignalParticle_charge", &MCSignalParticle_charge);
      output_tree->Branch("MCSignalParticle_Tauidx", &MCSignalParticle_Tauidx);
      output_tree->Branch("MCSignalParticle_SourceVertex", &MCSignalParticle_SourceVertex);
      output_tree->Branch("MCTauandProd_p4", &MCTauandProd_p4);
      output_tree->Branch("MCTauandProd_Vertex", &MCTauandProd_Vertex);
      output_tree->Branch("MCTauandProd_pdgid", &MCTauandProd_pdgid);
      output_tree->Branch("MCTauandProd_midx", &MCTauandProd_midx);
      output_tree->Branch("MCTauandProd_charge", &MCTauandProd_charge);
   }

   //================  Three Muonss block
   output_tree->Branch("ThreeMuons_index",&ThreeMuons_index);
   output_tree->Branch("ThreeMuons_SV_Chi2",&ThreeMuons_SV_Chi2);
   output_tree->Branch("ThreeMuons_SV_NDF",&ThreeMuons_SV_NDF);
   output_tree->Branch("ThreeMuons_TriggerMatch_dR",&ThreeMuons_TriggerMatch_dR);

   output_tree->Branch("signalTau_charge",&signalTau_charge);
   output_tree->Branch("signalTau_isLVP",&signalTau_isLVP);
   output_tree->Branch("signalTau_pdgid",&signalTau_pdgid);
   output_tree->Branch("signalTau_B", &signalTau_B);
   output_tree->Branch("signalTau_M",&signalTau_M);
   output_tree->Branch("signalTau_lvp",&signalTau_lvp);
   output_tree->Branch("signalTau_cov",&signalTau_cov);


   output_tree->Branch("TwoMuonsTrack_Muonsindex",&TwoMuonsTrack_Muonsindex);
   output_tree->Branch("TwoMuonsTrack_Trackindex",&TwoMuonsTrack_Trackindex);
   output_tree->Branch("TwoMuonsTrack_SV_Chi2",&TwoMuonsTrack_SV_Chi2);
   output_tree->Branch("TwoMuonsTrack_SV_NDF",&TwoMuonsTrack_SV_NDF);
   output_tree->Branch("TwoMuonsTrack_TriggerMatch_dR",&TwoMuonsTrack_TriggerMatch_dR);

   output_tree->Branch("Jet_BTagCVSB", &Jet_BTagCVSB);
   output_tree->Branch("Jet_BTagMVA", &Jet_BTagMVA);
   output_tree->Branch("Jet_BTagCSV", &Jet_BTagCSV);
   output_tree->Branch("Jet_p4",&Jet_p4);

   output_tree->Branch("Vertex_N_primary", &Vertex_N_primary);
   output_tree->Branch("Vertex_signal_dca_reco", &Vertex_signal_dca_reco);
   output_tree->Branch("Vertex_signal_KF_pos", &Vertex_signal_KF_pos);
   output_tree->Branch("Vertex_signal_KF_cov", &Vertex_signal_KF_cov);
   output_tree->Branch("Vertex_signal_KF_refittedTracksP4", &Vertex_signal_KF_refittedTracksP4);
   output_tree->Branch("Vertex_signal_KF_Chi2", &Vertex_signal_KF_Chi2);
   output_tree->Branch("Vertex_signal_AF_pos", &Vertex_signal_AF_pos);
   output_tree->Branch("Vertex_signal_AF_Chi2", &Vertex_signal_AF_Chi2);
   output_tree->Branch("Vertex_signal_AF_Ndf", &Vertex_signal_AF_Ndf);
   output_tree->Branch("Vertex_signal_KF_BS_2Ddistance", &Vertex_signal_KF_BS_2Ddistance);
   output_tree->Branch("Vertex_signal_KF_BS_error", &Vertex_signal_KF_BS_error);
   output_tree->Branch("Vertex_signal_KF_BS_significance", &Vertex_signal_KF_BS_significance);
   output_tree->Branch("Vertex_pair_quality", &Vertex_pair_quality);


   output_tree->Branch("Vertex_2MuonsIsoTrack_KF_Chi2", &Vertex_2MuonsIsoTrack_KF_Chi2 );
   output_tree->Branch("Vertex_2MuonsIsoTrack_KF_cov", &Vertex_2MuonsIsoTrack_KF_cov );
   output_tree->Branch("Vertex_2MuonsIsoTrack_KF_pos", &Vertex_2MuonsIsoTrack_KF_pos );



   output_tree->Branch("Vertex_Pair12_Pos", &Vertex_Pair12_Pos);
   output_tree->Branch("Vertex_Pair23_Pos", &Vertex_Pair23_Pos);
   output_tree->Branch("Vertex_Pair31_Pos", &Vertex_Pair31_Pos);

   output_tree->Branch("Vertex_pairfit_status", &Vertex_pairfit_status);
   output_tree->Branch("Vertex_MatchedPrimaryVertex",&Vertex_MatchedPrimaryVertex);
   output_tree->Branch("Vertex_SecondBestPrimaryVertex",&Vertex_SecondBestPrimaryVertex);
   output_tree->Branch("Vertex_RefitPVisValid",&Vertex_RefitPVisValid);
   output_tree->Branch("Vertex_MatchedRefitPrimaryVertex",&Vertex_MatchedRefitPrimaryVertex);
   output_tree->Branch("Vertex_MatchedRefitPrimaryVertex_covariance",&Vertex_MatchedRefitPrimaryVertex_covariance);
   output_tree->Branch("Vertex_HighestPt_PrimaryVertex",&Vertex_HighestPt_PrimaryVertex);
   output_tree->Branch("Vertex_HighestPt_PrimaryVertex_covariance",&Vertex_HighestPt_PrimaryVertex_covariance);
   output_tree->Branch("Vertex_d0_reco",&Vertex_d0_reco);
   output_tree->Branch("Vertex_dz_reco",&Vertex_dz_reco);
   output_tree->Branch("Vertex_d0SV_reco",&Vertex_d0SV_reco);
   output_tree->Branch("Vertex_dzSV_reco",&Vertex_dzSV_reco);
   output_tree->Branch("Vertex_d0BeamSpot_reco",&Vertex_d0BeamSpot_reco);
   output_tree->Branch("Vertex_d0BeamSpot_reco_sig",&Vertex_d0BeamSpot_reco_sig);
   output_tree->Branch("Vertex_d0sig_reco",&Vertex_d0sig_reco);
   output_tree->Branch("Vertex_d0sigSV_reco",&Vertex_d0sigSV_reco);
   output_tree->Branch("Vertex_2Ddisplacement",&Vertex_2Ddisplacement);
   output_tree->Branch("Vertex_3Ddisplacement",&Vertex_3Ddisplacement);
   output_tree->Branch("Vertex_Isolation1",&Vertex_Isolation1);
   output_tree->Branch("Vertex_Isolation2",&Vertex_Isolation2);
   output_tree->Branch("Vertex_Isolation3",&Vertex_Isolation3);
   output_tree->Branch("Vertex_Isolation4",&Vertex_Isolation4);

   output_tree->Branch("Vertex_NMuonsAssocWithPV",&Vertex_NMuonsAssocWithPV);


   output_tree->Branch("TriggerObject_pt",&TriggerObject_pt);
   output_tree->Branch("TriggerObject_phi",&TriggerObject_phi);
   output_tree->Branch("TriggerObject_eta",&TriggerObject_eta);
   output_tree->Branch("TriggerObject_name",&TriggerObject_name);

   //  output_tree->Branch("IsolationBranch_Trackp4", &IsolationBranch_Trackp4);



   output_tree->Branch("IsolationTrack_p4", &IsolationTrack_p4);

   output_tree->Branch("IsolationTrack_VertexWithSignalMuonIsValid", &IsolationTrack_VertexWithSignalMuonIsValid);
   output_tree->Branch("IsolationTrack_VertexWithSignalMuonChi2", &IsolationTrack_VertexWithSignalMuonChi2);
   output_tree->Branch("IsolationTrack_VertexWithSignalMuonPosition", &IsolationTrack_VertexWithSignalMuonPosition);


   output_tree->Branch("IsolationTrack_charge",&IsolationTrack_charge);
   output_tree->Branch("IsolationTrack_isHighPurity",&IsolationTrack_isHighPurity);
   output_tree->Branch("IsolationTrack_quality",&IsolationTrack_quality);
   output_tree->Branch("IsolationTrack_dxySV",&IsolationTrack_dxySV);
   output_tree->Branch("IsolationTrack_dzSV",&IsolationTrack_dzSV);
   output_tree->Branch("IsolationTrack_dxyPV",&IsolationTrack_dxyPV);
   output_tree->Branch("IsolationTrack_dzPV",&IsolationTrack_dzPV);
   output_tree->Branch("IsolationTrack_DocaMu1",&IsolationTrack_DocaMu1);
   output_tree->Branch("IsolationTrack_DocaMu2",&IsolationTrack_DocaMu2);
   output_tree->Branch("IsolationTrack_DocaMu3",&IsolationTrack_DocaMu3);


   output_tree->Branch("IsolationTrack_Helcharge", &IsolationTrack_Helcharge);
   output_tree->Branch("IsolationTrack_pdgid", &IsolationTrack_pdgid);
   output_tree->Branch("IsolationTrack_B", &IsolationTrack_B);
   output_tree->Branch("IsolationTrack_M", &IsolationTrack_M);
   output_tree->Branch("IsolationTrack_par", &IsolationTrack_par);
   output_tree->Branch("IsolationTrack_cov", &IsolationTrack_cov);




   output_tree->Branch("SV_Track_P4",&SV_Track_P4);
   output_tree->Branch("SV_pos",&SV_pos);
   output_tree->Branch("SV_TrackCharge",&SV_TrackCharge);
   output_tree->Branch("SV_Mass",&SV_Mass);
   output_tree->Branch("SV_PosCovariance",&SV_PosCovariance);



   output_tree->Branch("Trigger_l1name",&Trigger_l1name);
   output_tree->Branch("Trigger_l1decision",&Trigger_l1decision);
   output_tree->Branch("Trigger_l1prescale",&Trigger_l1prescale);

   output_tree->Branch("Trigger_hltname",&Trigger_hltname);
   output_tree->Branch("Trigger_hltdecision",&Trigger_hltdecision);

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



   dump_track_index_to_fill.clear();
   dump_pv_index_to_fill.clear();


   Event_nsignal_candidates=0;
   Event_ndsphipi_candidate=0;



   //=======  Taus  ===
   Tau_p4.clear();
   Tau_charge.clear();
   Tau_DecayMode.clear();
   Tau_DecayModeFinding.clear();
   Tau_NewDecayModeFinding.clear();


   Tau_byLooseDeepTau2017v2p1VSe.clear();
   Tau_byMediumDeepTau2017v2p1VSe.clear();
   Tau_byTightDeepTau2017v2p1VSe.clear();

   Tau_byLooseDeepTau2017v2p1VSmu.clear();
   Tau_byMediumDeepTau2017v2p1VSmu.clear();
   Tau_byTightDeepTau2017v2p1VSmu.clear();

   Tau_byLooseDeepTau2017v2p1VSjet.clear();
   Tau_byMediumDeepTau2017v2p1VSjet.clear();
   Tau_byTightDeepTau2017v2p1VSjet.clear();


   Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
   Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
   Tau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();
   Tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();

   Tau_FloatDiscriminants.clear();
   Tau_IntDiscriminants.clear();

   Tau_PFTauTrack_p4.clear();
   Tau_Track_par.clear();
   Tau_Track_cov.clear();

   Tau_Track_Charge.clear();
   Tau_Track_pdgid.clear();
   Tau_Track_B.clear();
   Tau_Track_M.clear();

   Tau_SVPos.clear();
   Tau_SVCov.clear();
   Tau_a1_charge.clear();
   Tau_a1_pdgid.clear();
   Tau_a1_B.clear();
   Tau_a1_M.clear();
   Tau_a1_lvp.clear();
   Tau_a1_cov.clear();




   //=======  Gammas ===
   Gamma_P4.clear();
   Gamma_hasPixelSeed.clear();
   Gamma_hasConversionTracks.clear();

   Gamma_e1x5.clear();
   Gamma_e2x5.clear();
   Gamma_e3x3.clear();
   Gamma_e5x5.clear();

   Gamma_isPFPhoton.clear();


   //=======  Electronss ===
   Electron_p4.clear();
   Electron_Charge.clear();
   Electron_puppiNeutralHadronIso.clear();
   Electron_puppiPhotonIso.clear();
   Electron_puppiChargedHadronIso.clear();
   Electron_trackIso.clear();
   Electron_isPF.clear();
   Electron_cutBasedElectronID_Fall17_94X_V2_veto.clear();
   Electron_cutBasedElectronID_Fall17_94X_V2_loose.clear();
   Electron_cutBasedElectronID_Fall17_94X_V2_medium.clear();
   Electron_cutBasedElectronID_Fall17_94X_V2_tight.clear();


   //=======  Muons ===
   Muon_p4.clear();
   Muon_Poca.clear();
   Muon_isGlobalMuon.clear();
   Muon_isPFMuon.clear();
   Muon_isRPCMuon.clear();
   Muon_isStandAloneMuon.clear();
   Muon_isTrackerMuon.clear();
   Muon_isCaloMuon.clear();
   Muon_isIsolationValid.clear();
   Muon_isQualityValid.clear();
   Muon_isTimeValid.clear();

   Muon_expectedNnumberOfMatchedStations.clear();


   Muon_emEt03.clear();
   Muon_emVetoEt03.clear();
   Muon_hadEt03.clear();
   Muon_hadVetoEt03.clear();
   Muon_nJets03.clear();
   Muon_nTracks03.clear();
   Muon_sumPt03.clear();
   Muon_trackerVetoPt03.clear();
   Muon_ID.clear();
   Muon_StandardSelection.clear();
   Muon_emEt05.clear();
   Muon_emVetoEt05.clear();
   Muon_hadEt05.clear();
   Muon_hadVetoEt05.clear();
   Muon_nJets05.clear();
   Muon_nTracks05.clear();
   Muon_sumPt05.clear();
   Muon_trackerVetoPt05.clear();
   Muon_timeAtIpInOut.clear();
   Muon_timeAtIpInOutErr.clear();


   Muon_sumChargedHadronPt03.clear();
   Muon_sumChargedParticlePt03.clear();
   Muon_sumNeutralHadronEt03.clear();
   Muon_sumNeutralHadronEtHighThreshold03.clear();
   Muon_sumPhotonEt03.clear();
   Muon_sumPhotonEtHighThreshold03.clear();
   Muon_sumPUPt03.clear();

   Muon_sumChargedHadronPt04.clear();
   Muon_sumChargedParticlePt04.clear();
   Muon_sumNeutralHadronEt04.clear();
   Muon_sumNeutralHadronEtHighThreshold04.clear();
   Muon_sumPhotonEt04.clear();
   Muon_sumPhotonEtHighThreshold04.clear();
   Muon_sumPUPt04.clear();

   Muon_numberOfChambers.clear();
   Muon_Track_idx.clear();

   Muon_charge.clear();
   Muon_trackCharge.clear();
   Muon_pdgid.clear();
   Muon_B.clear();
   Muon_M.clear();
   Muon_par.clear();
   Muon_cov.clear();



   Muon_hitPattern_pixelLayerwithMeas.clear();
   Muon_numberOfMatchedStations.clear();
   Muon_normChi2.clear();
   Muon_hitPattern_numberOfValidMuonHits.clear();
   Muon_innerTrack_numberofValidHits.clear();
   Muon_numberOfMatches.clear();
   Muon_numberofValidPixelHits.clear();
   Muon_trackerLayersWithMeasurement.clear();

   Muon_vmuonhitcomb_reco.clear();
   Muon_rpchits_reco.clear();

   Muon_combinedQuality_updatedSta.clear();
   Muon_combinedQuality_trkKink.clear();
   Muon_combinedQuality_glbKink.clear();
   Muon_combinedQuality_trkRelChi2.clear();
   Muon_combinedQuality_staRelChi2.clear();
   Muon_combinedQuality_chi2LocalPosition.clear();
   Muon_combinedQuality_chi2LocalMomentum.clear();
   Muon_combinedQuality_localDistance.clear();
   Muon_combinedQuality_globalDeltaEtaPhi.clear();
   Muon_combinedQuality_tightMatch.clear();
   Muon_combinedQuality_glbTrackProbability.clear();

   Muon_prod_inner_outer_charge.clear();
   Muon_outerTrack_p4.clear();
   Muon_innerTrack_p4.clear();
   Muon_innerTrack_quality.clear();
   Muon_ptErrOverPt.clear();
   Muon_calEnergy_hadS9.clear();
   Muon_calEnergy_had.clear();
   Muon_calEnergy_emS25.clear();
   Muon_calEnergy_emS9.clear();
   Muon_calEnergy_em.clear();

   Muon_segmentCompatibility.clear();
   Muon_caloCompatibility.clear();

   Muon_TrackX.clear();
   Muon_TrackY.clear();
   Muon_dDxDz.clear();
   Muon_dDyDz.clear();
   Muon_dX.clear();
   Muon_dY.clear();
   Muon_pullX.clear();
   Muon_pullY.clear();
   Muon_pullDxDz.clear();
   Muon_pullDyDz.clear();
   numberOfSegments.clear();


   Muon_ptError.clear();
   Muon_phiError.clear();
   Muon_etaError.clear();

   Muon_innerTrack_validFraction.clear();
   Muon_innerTrack_pixelLayersWithMeasurement.clear();
   Muon_innerTrack_numberOfValidTrackerHits.clear();
   Muon_innerTrack_numberOfLostTrackerHits.clear();
   Muon_innerTrack_numberOfLostTrackerInnerHits.clear();
   Muon_innerTrack_numberOfLostTrackerOuterHits.clear();
   Muon_innerTrack_normalizedChi2.clear();

   Muon_outerTrack_normalizedChi2.clear();
   Muon_outerTrack_muonStationsWithValidHits.clear();
   Muon_isGoodMuon_TM2DCompatibility.clear();
   Muon_isGoodMuon_TrackerMuonArbitrated.clear();
   Muon_isGoodMuon_TMOneStationTight.clear();
   Muon_isGoodMuon_TMOneStationAngTight.clear();
   Muon_isGoodMuon_TMLastStationTight.clear();
   Muon_isGoodMuon_TMLastStationAngTight.clear();
   Muon_isGoodMuon_TMLastStationOptimizedLowPtTight.clear();
   Muon_isGoodMuon_TMLastStationOptimizedBarrelLowPtTight.clear();

   TriggerObject_pt.clear();
   TriggerObject_eta.clear();
   TriggerObject_phi.clear();
   TriggerObject_name.clear();

   if (doMC_) {
      MC_isReco=0;
      MC_p4.clear();
      MC_vertex.clear();
      MC_pdgid.clear();
      MC_charge.clear();
      MC_midx.clear();
      MC_status.clear();
      MC_childpdgid.clear();
      MC_childidx.clear();
      MCSignalParticle_p4.clear();
      MCSignalParticle_Vertex.clear();
      MCSignalParticle_pdgid.clear();
      MCSignalParticle_charge.clear();
      MCSignalParticle_Poca.clear();
      MCSignalParticle_Tauidx.clear();
      MCTauandProd_p4.clear();
      MCTauandProd_Vertex.clear();
      MCTauandProd_pdgid.clear();
      MCTauandProd_midx.clear();
      MCTauandProd_charge.clear();
      MCTau_JAK.clear();
      MCTau_DecayBitMask.clear();
      MCSignalParticle_childpdgid.clear();
      MCSignalParticle_childp4.clear();
      MCSignalParticle_Sourcepdgid.clear();
      MCSignalParticle_Sourcep4.clear();
      MCSignalParticle_SourceVertex.clear();

   }


   ThreeMuons_idx.clear();
   ThreeMuons_index.clear();
   ThreeMuons_SV_Chi2.clear();
   ThreeMuons_SV_NDF.clear();
   ThreeMuons_TriggerMatch_dR.clear();

   signalTau_charge.clear();
   signalTau_isLVP.clear();
   signalTau_pdgid.clear();
   signalTau_B.clear();
   signalTau_M.clear();
   signalTau_lvp.clear();
   signalTau_cov.clear();

   TwoMuonsTrack_idx.clear();
   TwoMuonsTrack_Muonsindex.clear();
   TwoMuonsTrack_Trackindex.clear();
   TwoMuonsTrack_SV_Chi2.clear();
   TwoMuonsTrack_SV_NDF.clear();
   TwoMuonsTrack_TriggerMatch_dR.clear();

   Jet_BTagCVSB.clear();
   Jet_BTagMVA.clear();
   Jet_BTagCSV.clear();
   Jet_p4.clear();

   Vertex_NMuonsAssocWithPV.clear();
   Vertex_signal_dca_reco.clear();
   Vertex_signal_KF_pos.clear();
   Vertex_signal_KF_cov.clear();
   Vertex_signal_KF_refittedTracksP4.clear();
   Vertex_signal_KF_Chi2.clear();
   Vertex_signal_AF_pos.clear();
   Vertex_signal_AF_Chi2.clear();
   Vertex_signal_AF_Ndf.clear();
   Vertex_signal_KF_BS_2Ddistance.clear();
   Vertex_signal_KF_BS_error.clear();
   Vertex_signal_KF_BS_significance.clear();

   Vertex_2MuonsIsoTrack_KF_Chi2.clear();
   Vertex_2MuonsIsoTrack_KF_cov.clear();
   Vertex_2MuonsIsoTrack_KF_pos.clear();



   Vertex_pair_quality.clear();
   Vertex_pairfit_status.clear();
   Vertex_Pair12_Pos.clear();
   Vertex_Pair23_Pos.clear();
   Vertex_Pair31_Pos.clear();

   Vertex_MatchedPrimaryVertex.clear();
   Vertex_SecondBestPrimaryVertex.clear();

   Vertex_RefitPVisValid.clear();
   Vertex_MatchedRefitPrimaryVertex.clear();
   Vertex_MatchedRefitPrimaryVertex_covariance.clear();

   Vertex_HighestPt_PrimaryVertex.clear();
   Vertex_HighestPt_PrimaryVertex_covariance.clear();


   Vertex_d0_reco.clear();
   Vertex_dz_reco.clear();
   Vertex_d0SV_reco.clear();
   Vertex_dzSV_reco.clear();
   Vertex_d0sig_reco.clear();
   Vertex_d0sigSV_reco.clear();
   Vertex_d0BeamSpot_reco.clear();
   Vertex_d0BeamSpot_reco_sig.clear();
   Vertex_2Ddisplacement.clear();
   Vertex_3Ddisplacement.clear();
   Vertex_Isolation1.clear();
   Vertex_Isolation2.clear();
   Vertex_Isolation3.clear();
   Vertex_Isolation4.clear();

   SV_Track_P4.clear();
   SV_pos.clear();
   SV_Mass.clear();
   SV_PosCovariance.clear();
   SV_TrackCharge.clear();

   IsolationBranch_Trackp4.clear();

   IsolationTrack_p4.clear();
   IsolationTrack_VertexWithSignalMuonIsValid.clear();
   IsolationTrack_VertexWithSignalMuonChi2.clear();
   IsolationTrack_VertexWithSignalMuonPosition.clear();
   IsolationTrack_charge.clear();
   IsolationTrack_isHighPurity.clear();
   IsolationTrack_quality.clear();
   IsolationTrack_dxySV.clear();
   IsolationTrack_dzSV.clear();
   IsolationTrack_dxyPV.clear();
   IsolationTrack_dzPV.clear();
   IsolationTrack_DocaMu1.clear();
   IsolationTrack_DocaMu2.clear();
   IsolationTrack_DocaMu3.clear();


   IsolationTrack_Helcharge.clear();
   IsolationTrack_pdgid.clear();
   IsolationTrack_B.clear();
   IsolationTrack_M.clear();
   IsolationTrack_par.clear();
   IsolationTrack_cov.clear();



   Trigger_l1name.clear();
   Trigger_l1decision.clear();
   Trigger_l1prescale.clear();
   Trigger_hltname.clear();
   Trigger_hltdecision.clear();
   
   tauIntDiscrims_.clear();
   tauFloatDiscrims_.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(T3MNtuple);

