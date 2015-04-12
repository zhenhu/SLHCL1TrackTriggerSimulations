#include "SLHCL1TrackTriggerSimulations/NTupleTools/interface/NTupleGenParticlesExtra.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


NTupleGenParticlesExtra::NTupleGenParticlesExtra(const edm::ParameterSet& iConfig) :
  inputTag_(iConfig.getParameter<edm::InputTag>("inputTag")),
  prefix_  (iConfig.getParameter<std::string>("prefix")),
  suffix_  (iConfig.getParameter<std::string>("suffix")),
  selector_(iConfig.existsAs<std::string>("cut") ? iConfig.getParameter<std::string>("cut") : "", true),
  maxN_    (iConfig.getParameter<unsigned>("maxN")) {

    //produces<std::vector<float> > (prefix_ + "phi"      + suffix_);

    produces<std::vector<float> > (prefix_ + "qoverp"   + suffix_);
    produces<std::vector<float> > (prefix_ + "lambda"   + suffix_);
    produces<std::vector<float> > (prefix_ + "dxy"      + suffix_);
    produces<std::vector<float> > (prefix_ + "dsz"      + suffix_);

    produces<std::vector<float> > (prefix_ + "invPt"    + suffix_);
    produces<std::vector<float> > (prefix_ + "cotTheta" + suffix_);
    produces<std::vector<float> > (prefix_ + "d0"       + suffix_);
    produces<std::vector<float> > (prefix_ + "dz"       + suffix_);
}

void NTupleGenParticlesExtra::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    //std::auto_ptr<std::vector<float> > v_phi      (new std::vector<float>());

    std::auto_ptr<std::vector<float> > v_qoverp   (new std::vector<float>());
    std::auto_ptr<std::vector<float> > v_lambda   (new std::vector<float>());
    std::auto_ptr<std::vector<float> > v_dxy      (new std::vector<float>());
    std::auto_ptr<std::vector<float> > v_dsz      (new std::vector<float>());

    std::auto_ptr<std::vector<float> > v_invPt    (new std::vector<float>());
    std::auto_ptr<std::vector<float> > v_cotTheta (new std::vector<float>());
    std::auto_ptr<std::vector<float> > v_d0       (new std::vector<float>());
    std::auto_ptr<std::vector<float> > v_dz       (new std::vector<float>());

    //__________________________________________________________________________
    if (!iEvent.isRealData()) {
        edm::Handle<reco::GenParticleCollection> parts;
        iEvent.getByLabel(inputTag_, parts);

        if (parts.isValid()) {
            edm::LogInfo("NTupleGenParticlesExtra") << "Size: " << parts->size();

            unsigned n = 0;
            for (reco::GenParticleCollection::const_iterator it = parts->begin(); it != parts->end(); ++it) {
                if (n >= maxN_)
                    break;
                if (!selector_(*it))
                    continue;

                double theta    = it->theta();
                double lambda   = M_PI_2 - theta;
                double cotTheta = 1.0 / std::tan(theta);

                // Fill the vectors
                //v_phi      ->push_back(it->phi());

                v_qoverp   ->push_back(float(it->charge()) / it->p());
                v_lambda   ->push_back(lambda);
                //v_dxy      ->push_back(dxy);
                //v_dsz      ->push_back(dsz);

                v_invPt    ->push_back(float(it->charge()) / it->pt());
                v_cotTheta ->push_back(cotTheta);
                //v_d0       ->push_back(d0);
                //v_dz       ->push_back(dz);

                n++;
            }

        } else {
            edm::LogError("NTupleGenParticlesExtra") << "Cannot get the product: " << inputTag_;
        }
    }

    //__________________________________________________________________________
    //iEvent.put(v_phi      , prefix_ + "phi"      + suffix_);

    iEvent.put(v_qoverp   , prefix_ + "qoverp"   + suffix_);
    iEvent.put(v_lambda   , prefix_ + "lambda"   + suffix_);
    iEvent.put(v_dxy      , prefix_ + "dxy"      + suffix_);
    iEvent.put(v_dsz      , prefix_ + "dsz"      + suffix_);

    iEvent.put(v_invPt    , prefix_ + "invPt"    + suffix_);
    iEvent.put(v_cotTheta , prefix_ + "cotTheta" + suffix_);
    iEvent.put(v_d0       , prefix_ + "d0"       + suffix_);
    iEvent.put(v_dz       , prefix_ + "dz"       + suffix_);
}
