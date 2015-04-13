#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoPCA.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulation/external/Eigen/Eigenvalues"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/external/Eigen/QR"
#include <iostream>
#include <iomanip>
#include <fstream>

static const unsigned NVARIABLES = 12;  // number of hit coordinates
static const unsigned NPARAMETERS = 4;  // number of track parameters


// _____________________________________________________________________________
int TrackFitterAlgoPCA::loadConstants(TString txt) {
    std::ifstream infile(txt.Data());
    if (!infile) {
        std::cout << "Unable to open " << txt << std::endl;
        return 1;
    }

    double x;
    V_ = Eigen::MatrixXd::Zero(NVARIABLES, NVARIABLES);
    for (unsigned ivar=0; ivar<NVARIABLES; ++ivar) {
        for (unsigned jvar=0; jvar<NVARIABLES; ++jvar) {
            infile >> x;
            V_(ivar, jvar) = x;
        }
    }

    D_ = Eigen::MatrixXd::Zero(NVARIABLES, NPARAMETERS);
    for (unsigned ipar=0; ipar<NPARAMETERS; ++ipar) {
        for (unsigned jvar=0; jvar<NVARIABLES; ++jvar) {
            infile >> x;
            D_(ipar, jvar) = x;
        }
    }

    DV_ = Eigen::MatrixXd::Zero(NPARAMETERS, NVARIABLES);
    DV_ = D_ * V_;

    return 0;
}

// _____________________________________________________________________________
int TrackFitterAlgoPCA::fit(const std::vector<TTHit>& hits, TTTrack2& track) {
    return 0;
}

// _____________________________________________________________________________
void TrackFitterAlgoPCA::print() {
    std::ios::fmtflags flags = std::cout.flags();
    std::cout << std::setprecision(4);
    std::cout << "V: " << std::endl;
    std::cout << V_ << std::endl << std::endl;
    std::cout << "D: " << std::endl;
    std::cout << D_ << std::endl << std::endl;
    std::cout << "DV: " << std::endl;
    std::cout << DV_ << std::endl << std::endl;
    std::cout.flags(flags);
}