#ifndef AMSimulation_TrackFitterAlgoPCA_h_
#define AMSimulation_TrackFitterAlgoPCA_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoBase.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PCA.h"
using namespace slhcl1tt;

#include "SLHCL1TrackTriggerSimulations/AMSimulation/external/Eigen/Core"


class TrackFitterAlgoPCA : public TrackFitterAlgoBase {
  public:
    TrackFitterAlgoPCA(const slhcl1tt::ProgramOption& po)
    : TrackFitterAlgoBase() {

        if (po.view == "XYZ" || po.view == "3D")
            view_ = PCA_3D;
        else if (po.view == "XY" || po.view == "RPHI")
            view_ = PCA_RPHI;
        else if (po.view == "RZ")
            view_ = PCA_RZ;

        if (po.algo == "PCA4")
            fiveParams_ = false;
        else
            fiveParams_ = true;

        hitbits_ = static_cast<PCA_HitBits>(po.hitbits);

        // Book histograms
        bookHistograms();
    }

    ~TrackFitterAlgoPCA() {}

    int bookHistograms();

    int loadConstants(TString txt);

    int fit(const std::vector<TTHit>& hits, TTTrack2& track);

    void print();

  private:
    // Settings
    bool fiveParams_;
    PCA_FitView view_;
    PCA_HitBits hitbits_;

    // Matrices
    Eigen::VectorXd shifts_;
    Eigen::VectorXd sqrtEigenvalues_;
    Eigen::MatrixXd D_;
    Eigen::MatrixXd V_;
    Eigen::MatrixXd DV_;
};

#endif
