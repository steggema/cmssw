#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include <vector>
#include <iostream>


class PFCandidateJetLinker : public edm::global::EDProducer<> {
    public:
        PFCandidateJetLinker( edm::ParameterSet const & params ) :
            objName_(params.getParameter<std::string>("objName")),
            name_(params.getParameter<std::string>("name")),
            doc_(params.getParameter<std::string>("doc")),
            pfcandidates_(consumes<edm::View<reco::Candidate>>(params.getParameter<edm::InputTag>("pfcandidates"))),
            jets_(consumes<edm::View<pat::Jet>>(params.getParameter<edm::InputTag>("jets")))
        {
            produces<nanoaod::FlatTable>(name_); // "PF"
        }

        ~PFCandidateJetLinker() override {}

        void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override {

            edm::Handle<edm::View<pat::Jet>> jets;
            iEvent.getByToken(jets_, jets);

            edm::Handle<edm::View<reco::Candidate>> pfcandidates;
            iEvent.getByToken(pfcandidates_, pfcandidates);
            unsigned int ncand = pfcandidates->size();
            std::vector<int> pfcand2jet_index(ncand, -1);

            // loop over pfcandidates, then match with jets
            int jet_index = 0;
            for (const auto& jet : *jets) {
                auto constituents = jet.daughterPtrVector();
                for (const auto& jetconst : constituents) {
                        pfcand2jet_index[jetconst.key()] = jet_index;
                }
                ++jet_index;
            }

            auto tab = std::make_unique<nanoaod::FlatTable>(ncand, name_, false, true);

            tab->addColumn<int>(objName_ + "Idx_cool", pfcand2jet_index, "Index into " + objName_ + " list for " + name_, nanoaod::FlatTable::IntColumn);

            iEvent.put(std::move(tab), name_);
        }


        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
            edm::ParameterSetDescription desc;
            desc.add<std::string>("objName")->setComment("name of the nanoaod::FlatTable to extend with this table");
            desc.add<std::string>("name")->setComment("name of the nanoaod::FlatTable linked to");
            desc.add<std::string>("doc")->setComment("Add Jet indices to PFCandidate collection");
            desc.add<edm::InputTag>("pfcandidates")->setComment("pfcandidates InputTag");
            desc.add<edm::InputTag>("jets")->setComment("jets InputTag");
            descriptions.add("pfCandidateJetLinker", desc);
        }

    protected:
        const std::string objName_, name_, doc_;
        const edm::EDGetTokenT<edm::View<reco::Candidate>> pfcandidates_;
        const edm::EDGetTokenT<edm::View<pat::Jet>> jets_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFCandidateJetLinker);

