#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"

void T3MNtuple::fillMCTruth(const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      const Handle<GenParticleCollection>& genParticles,
      const Handle<vector<PileupSummaryInfo>>& puInfo)
{
  //  std::cout<<"do we fill MC  ??? "<< std::endl;
  if (!iEvent.isRealData())
    {
      if(doFullMC_){

         for (reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr != genParticles->end(); ++itr) {
            if(SkipThisParticle(*itr)) continue;
            MC_pdgid.push_back(itr->pdgId());
            MC_charge.push_back(itr->charge());
            std::vector<float> iMC_p4;
            std::vector<float> iMC_vertex;
            iMC_p4.push_back(itr->p4().E());
            iMC_p4.push_back(itr->p4().Px());
            iMC_p4.push_back(itr->p4().Py());
            iMC_p4.push_back(itr->p4().Pz());

            iMC_vertex.push_back(itr->vx());
            iMC_vertex.push_back(itr->vy());
            iMC_vertex.push_back(itr->vz());

            MC_p4.push_back(iMC_p4);
            MC_vertex.push_back(iMC_vertex);
            MC_midx.push_back(-1);
            MC_status.push_back(itr->status());
            MC_childpdgid.push_back(std::vector<int>());
            MC_childidx.push_back(std::vector<int>());
         }
         unsigned int i = 0;
         for (reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr != genParticles->end(); ++itr) {
            if(SkipThisParticle(*itr)) continue;
            for (unsigned int d = 0; d < itr->numberOfDaughters(); d++) {
               const reco::GenParticle *dau = static_cast<const reco::GenParticle*>(itr->daughter(d));
               unsigned int j = 0;
               for (reco::GenParticleCollection::const_iterator jtr = genParticles->begin(); jtr != genParticles->end(); ++jtr){
                  if(SkipThisParticle(*jtr)) continue;
                  if (dau->status() == jtr->status() && dau->p4() == jtr->p4() && dau->pdgId() == jtr->pdgId() && dau->numberOfMothers() == jtr->numberOfMothers()
                        && dau->numberOfDaughters() == jtr->numberOfDaughters()) {
                     MC_midx.at(j) = i;
                     MC_childidx.at(i).push_back(j);
                     MC_childpdgid.at(i).push_back(dau->pdgId());
                  }
                  j++;
               }
            }
            i++;
         }
      }

      TauDecay_CMSSW myTauDecay;
      DataMCType DMT;
      DataMC_Type_idx = DMT.GetType();
      myTauDecay.CheckForSignal(DataMC_Type_idx, genParticles);
      Event_DataMC_Type=DataMC_Type_idx;

      //genParticles
      unsigned int k(0);
      for (reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr != genParticles->end(); ++itr) {
         if (DMT.isSignalParticle(itr->pdgId())) {
            MCSignalParticle_childpdgid.push_back(std::vector<int>());
            MCSignalParticle_childp4.push_back(std::vector<std::vector<float> >());
            MCSignalParticle_Sourcepdgid.push_back(std::vector<int>());
            MCSignalParticle_Sourcep4.push_back(std::vector<std::vector<float> >());
            MCSignalParticle_pdgid.push_back(itr->pdgId());
            MCSignalParticle_charge.push_back(itr->charge());
            MCSignalParticle_Tauidx.push_back(std::vector<unsigned int>());

            std::vector<float> iSig_p4;
            iSig_p4.push_back(itr->p4().E());
            iSig_p4.push_back(itr->p4().Px());
            iSig_p4.push_back(itr->p4().Py());
            iSig_p4.push_back(itr->p4().Pz());
            MCSignalParticle_p4.push_back(iSig_p4);

            std::vector<float> iSig_Vertex;
            iSig_Vertex.push_back(itr->vx());
            iSig_Vertex.push_back(itr->vy());
            iSig_Vertex.push_back(itr->vz());
            MCSignalParticle_Vertex.push_back(iSig_Vertex);

            std::vector<float> iSourceVtx;

            if(itr->numberOfMothers()!=0){
               iSourceVtx.push_back(itr->mother(0)->vx());
               iSourceVtx.push_back(itr->mother(0)->vy());
               iSourceVtx.push_back(itr->mother(0)->vz());
            }
            MCSignalParticle_SourceVertex.push_back(iSourceVtx);

            for (unsigned int i = 0; i < itr->numberOfMothers(); i++){
               const reco::Candidate *mot = itr->mother(i);
               std::vector<float> iSourcep4;
               iSourcep4.push_back(mot->p4().E());
               iSourcep4.push_back(mot->p4().Px());
               iSourcep4.push_back(mot->p4().Py());
               iSourcep4.push_back(mot->p4().Pz());



               MCSignalParticle_Sourcepdgid.at(MCSignalParticle_Sourcepdgid.size() - 1).push_back(mot->pdgId());
               MCSignalParticle_Sourcep4.at(MCSignalParticle_Sourcepdgid.size() - 1).push_back(iSourcep4);
            }

            // look for daughter tau
            for (unsigned int i = 0; i < itr->numberOfDaughters(); i++){
               const reco::Candidate *dau = itr->daughter(i);
               std::vector<float> ichildp4;
               ichildp4.push_back(dau->p4().E());
               ichildp4.push_back(dau->p4().Px());
               ichildp4.push_back(dau->p4().Py());
               ichildp4.push_back(dau->p4().Pz());

               MCSignalParticle_childpdgid.at(MCSignalParticle_childpdgid.size() - 1).push_back(dau->pdgId());
               MCSignalParticle_childp4.at(MCSignalParticle_childpdgid.size() - 1).push_back(ichildp4);
               if (abs(dau->pdgId()) == PDGInfo::tau_minus ) {
                  unsigned int tauidx = MCTauandProd_p4.size();
                  MCSignalParticle_Tauidx.at(MCSignalParticle_Tauidx.size() - 1).push_back(tauidx);
                  // Analysis the tau decay

		  unsigned int JAK_ID, TauBitMask;
		  myTauDecay.AnalyzeTau(static_cast<const reco::GenParticle*>(dau), JAK_ID, TauBitMask);
                  std::vector<const reco::GenParticle*> TauProducts = TauDecayProducts(static_cast<const reco::GenParticle*>(dau));

                  MCTauandProd_midx.push_back(k);
                  MCTauandProd_pdgid.push_back(std::vector<int>());
                  MCTauandProd_charge.push_back(std::vector<int>());
                  MCTauandProd_p4.push_back(std::vector<std::vector<float> >());
                  MCTauandProd_Vertex.push_back(std::vector<std::vector<float> >());
                  for (unsigned int j = 0; j < TauProducts.size(); j++) {
		     MCTauandProd_pdgid.at(tauidx).push_back(TauProducts.at(j)->pdgId());
                     MCTauandProd_charge.at(tauidx).push_back(TauProducts.at(j)->charge());
                     std::vector<float> iTauandProd_p4;
                     std::vector<float> iTauandProd_vertex;
                     iTauandProd_p4.push_back(TauProducts.at(j)->p4().E());
                     iTauandProd_p4.push_back(TauProducts.at(j)->p4().Px());
                     iTauandProd_p4.push_back(TauProducts.at(j)->p4().Py());
                     iTauandProd_p4.push_back(TauProducts.at(j)->p4().Pz());

                     iTauandProd_vertex.push_back(TauProducts.at(j)->vx());
                     iTauandProd_vertex.push_back(TauProducts.at(j)->vy());
                     iTauandProd_vertex.push_back(TauProducts.at(j)->vz());
                     MCTauandProd_p4.at(tauidx).push_back(iTauandProd_p4);
                     MCTauandProd_Vertex.at(tauidx).push_back(iTauandProd_vertex);
                  }
               }
            }
            k++;
         }
      }

      puN = puInfo->begin()->getTrueNumInteractions();
   }
}
