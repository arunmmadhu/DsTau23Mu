#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"

void T3MNtuple::fillBTagJets( const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<JetTagCollection>& btagsCvsB,
      const Handle<JetTagCollection>& btagsCSV, 
      const Handle<JetTagCollection>& btagsMVA)
{

   for(size_t j = 0 ; j < btagsCvsB->size(); j++) {
      const JetTag & btag1 = (*btagsCvsB)[j];
      if(btag1.first->pt() > 5){
         std::vector<double> iJet_p4;
         iJet_p4.push_back(btag1.first->p4().e());
         iJet_p4.push_back(btag1.first->p4().px());
         iJet_p4.push_back(btag1.first->p4().py());
         iJet_p4.push_back(btag1.first->p4().pz());
         Jet_p4.push_back(iJet_p4);
         Jet_BTagCVSB.push_back(btag1.second);
         const JetTag & btag2 = (*btagsMVA)[j];
         Jet_BTagMVA.push_back(btag2.second);
         const JetTag & btag3 = (*btagsCSV)[j];
         Jet_BTagCSV.push_back(btag3.second<0 ? 0:btag3.second);
      }
   }
}
