#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"

void T3MNtuple::fillTrigger(const edm::Event& iEvent,
                            const edm::EventSetup& iSetup,
                            const Handle<TriggerResults>& triggerBitsH,
                            const Handle<trigger::TriggerEvent>& triggerSummary,
                            const Handle<vector<pat::TriggerObjectStandAlone> >& triggerObjects,
                            const TriggerNames& triggerNames)
{
   gtUtil_->retrieveL1(iEvent, iSetup, algToken_);
   const vector<pair<string, bool> > initialDecisions = gtUtil_->decisionsInitial();

   if (!iEvent.isRealData())
   {
      for (size_t i_l1t = 0; i_l1t < initialDecisions.size(); i_l1t++){

         string l1tName = (initialDecisions.at(i_l1t)).first;
         if(l1tName.find("DoubleMu") != string::npos || l1tName.find("TripleMu") != string::npos || l1tName.find("singleMu"))
         {
            Trigger_l1name.push_back( l1tName );

            Trigger_l1decision.push_back( initialDecisions.at(i_l1t).second );
            Trigger_l1prescale.push_back( 1 );
         }
      }
   }
   else
   {
      ESHandle<L1TGlobalPrescalesVetos> psAndVetos;
      auto psRcd = iSetup.tryToGet<L1TGlobalPrescalesVetosRcd>();
      if(psRcd) psRcd->get(psAndVetos);
      int columnN= gtUtil_->prescaleColumn();
      for (size_t i_l1t = 0; i_l1t < initialDecisions.size(); i_l1t++) {
         string l1tName = (initialDecisions.at(i_l1t)).first;
         if(l1tName.find("DoubleMu") != string::npos || l1tName.find("TripleMu") != string::npos || l1tName.find("SingleMu"))
         {
            Trigger_l1name.push_back( l1tName );
            Trigger_l1decision.push_back( initialDecisions.at(i_l1t).second );
            Trigger_l1prescale.push_back( (psAndVetos->prescale_table_)[columnN][i_l1t]);
         }
      }
   } 

   for (size_t i_hlt = 0; i_hlt != triggerBitsH->size(); ++i_hlt)
   {
      string hltName = triggerNames.triggerName(i_hlt);
      if(hltName.find("HLT_DoubleMu") != string::npos  or hltName.find("HLT_Mu") != string::npos or hltName.find("HLT_Dimuon0") != string::npos)
      {
         Trigger_hltname.push_back(hltName);
         Trigger_hltdecision.push_back(triggerBitsH->accept(i_hlt ));
      }
   }

   // Fill Trigger object info
   if (miniAODRun_){

      for (pat::TriggerObjectStandAlone to: *triggerObjects){

         to.unpackPathNames(triggerNames);
         trigger::size_type nFilters = to.filterLabels().size();
         for (trigger::size_type iFilter=0; iFilter<nFilters; iFilter++){
            std::string filterName = to.filterLabels()[iFilter];
            
            if (filterName.compare("hltTau3muTkVertexFilter")!=0 &&
                filterName.compare("hltdstau3muDisplaced3muFltr")!=0 &&
                //filterName.compare("hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p07")!=0 &&
                filterName.compare("hltDisplacedmumuFilterDimuon0LowMass")) continue;
            
            TriggerObject_name.push_back(filterName);
            TriggerObject_pt.push_back(to.pt());
            TriggerObject_eta.push_back(to.eta());
            TriggerObject_phi.push_back(to.phi());
         }
      }

   } else {

      std::vector<trigger::TriggerObject> trgobjs = triggerSummary->getObjects();

      edm::InputTag MuonFilterTag2017("hltTau3muTkVertexFilter","","HLT");
      edm::InputTag MuonFilterTag2018("hltdstau3muDisplaced3muFltr","","HLT");

      size_t MuonFilterIndex = (*triggerSummary).filterIndex(MuonFilterTag2017);
      if(MuonFilterIndex < (*triggerSummary).sizeFilters()) {
         const trigger::Keys &KEYS = (*triggerSummary).filterKeys(MuonFilterIndex);
         for (unsigned int ipart = 0; ipart < KEYS.size(); ipart++) {
            TriggerObject_name.push_back("hltTau3muTkVertexFilter");
            TriggerObject_pt.push_back(trgobjs.at(KEYS.at(ipart)).pt());
            TriggerObject_eta.push_back(trgobjs.at(KEYS.at(ipart)).eta());
            TriggerObject_phi.push_back(trgobjs.at(KEYS.at(ipart)).phi());
         }  
      }

      MuonFilterIndex = (*triggerSummary).filterIndex(MuonFilterTag2018);
      if(MuonFilterIndex < (*triggerSummary).sizeFilters()) {
         const trigger::Keys &KEYS = (*triggerSummary).filterKeys(MuonFilterIndex);
         for (unsigned int ipart = 0; ipart < KEYS.size(); ipart++) {
            TriggerObject_name.push_back("hltdstau3muDisplaced3muFltr");
            TriggerObject_pt.push_back(trgobjs.at(KEYS.at(ipart)).pt());
            TriggerObject_eta.push_back(trgobjs.at(KEYS.at(ipart)).eta());
            TriggerObject_phi.push_back(trgobjs.at(KEYS.at(ipart)).phi());
         }  
      }
   }
}
