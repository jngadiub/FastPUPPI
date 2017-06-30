#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "FastPUPPI/NtupleProducer/interface/L1TPFParticle.h"
#include "DataFormats/Math/interface/deltaPhi.h"

namespace l1tpf {
    class EcalProducerFromTPDigi : public edm::stream::EDProducer<> {
        public:
            explicit EcalProducerFromTPDigi(const edm::ParameterSet&) ;
            ~EcalProducerFromTPDigi() {}

        private:
            edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> EcalTPTag_;
            double etCut_;
            edm::ESHandle<CaloGeometry> pG;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;

            struct SimpleHit { 
	      float et, eta, phi; DetId id;
	      SimpleHit(float aet, float aeta, float aphi, DetId aid) : et(aet), eta(aeta), phi(aphi), id(aid) {}
            };

    }; // class
} // namespace

l1tpf::EcalProducerFromTPDigi::EcalProducerFromTPDigi(const edm::ParameterSet & iConfig) :
    EcalTPTag_(consumes<EcalEBTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("EcalTPTag"))),
    etCut_(iConfig.getParameter<double>("etMin"))

{
    produces<std::vector<l1tpf::Particle>>("crystals");
    produces<std::vector<l1tpf::Particle>>("towers");
}


void 
l1tpf::EcalProducerFromTPDigi::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
    std::unique_ptr<std::vector<l1tpf::Particle>> out_crystal(new std::vector<l1tpf::Particle>());
    std::unique_ptr<std::vector<l1tpf::Particle>> out_tower(new std::vector<l1tpf::Particle>());

    std::map<std::pair<int,int>,std::vector<SimpleHit>> towers;

    //Get Calo Geometry
    iSetup.get<CaloGeometryRecord>().get(pG);
    const CaloGeometry* caloGeom = pG.product();

    edm::Handle< EcalEBTrigPrimDigiCollection > ecalTPs;
    iEvent.getByToken(EcalTPTag_, ecalTPs);
    for (const EcalEBTriggerPrimitiveDigi & digi: *ecalTPs) {
        //https://github.com/cms-sw/cmssw/blob/master/SimCalorimetry/EcalEBTrigPrimProducers/plugins/EcalEBTrigPrimAnalyzer.cc#L173-L207
        const EBDetId TPid = digi.id();
        float et = digi.encodedEt()/8.; // 8 ADCcounts/GeV
        if (et <= etCut_) continue;
        const GlobalPoint & pos = caloGeom->getPosition(TPid);
        out_crystal->emplace_back(et, pos.eta(), pos.phi(), 0., 0);
	DetIdCollection detIds_;
	detIds_.push_back(TPid);
	out_crystal->back().setDetIds(detIds_);
	out_crystal->back().setSeed(TPid);
	towers[std::make_pair(TPid.tower_ieta(),TPid.tower_iphi())].emplace_back(et, pos.eta(), pos.phi(), TPid);
    }

    for (const auto & pair : towers) {
        double etsum = 0., etaetsum = 0., phietsum = 0.; DetIdCollection idsum;
        double reta = pair.second.front().eta, rphi = pair.second.front().phi;
        for (const SimpleHit & hit : pair.second) {
	  etsum += hit.et;
	  etaetsum += (hit.eta - reta) * hit.et;
	  phietsum += reco::deltaPhi(hit.phi, rphi) * hit.et;
	  idsum.push_back(hit.id);
        }
	if(etsum > 0)  {
	  etaetsum /= etsum;
	  phietsum /= etsum;
	}
	out_tower->emplace_back(etsum, etaetsum + reta, reco::deltaPhi(phietsum + rphi, 0.), 0., 0);
	out_tower->back().setDetIds(idsum);
	out_tower->back().setSeed(DetId(0));
    }

    iEvent.put(std::move(out_crystal), "crystals");
    iEvent.put(std::move(out_tower), "towers");
}
using l1tpf::EcalProducerFromTPDigi;
DEFINE_FWK_MODULE(EcalProducerFromTPDigi);
