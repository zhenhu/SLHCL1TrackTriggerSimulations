#ifndef AMSimulation_TrackFitterAlgoLTF_h_
#define AMSimulation_TrackFitterAlgoLTF_h_

#include <memory>
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoBase.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitter.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/LinearizedTrackFitterEmulator.h"

namespace slhcl1tt {

// The second parameter in the LinearizedTrackFitter constructor specifies if the fitter will replace the R coordinate of stubs in 2S modules in the disks
// with the extrapolated R using a tan(theta) pre-estimate.
class TrackFitterAlgoLTF : public TrackFitterAlgoBase
{
 public:
  TrackFitterAlgoLTF(const slhcl1tt::ProgramOption& po) :
    TrackFitterAlgoBase(),
    verbose_(po.verbose), emu_(po.emu), normChi2_(-1.)
  {
    if (po.emu == 0) linearizedTrackFitter_ = std::make_shared<LinearizedTrackFitter>("LinearizedTrackFit/LinearizedTrackFit/python/", true, 0, true, 14);
    else linearizedTrackFitterEmulator_ = std::make_shared<LinearizedTrackFitterEmulator>("LinearizedTrackFit/LinearizedTrackFit/python/", true, 0, true, 14);
  }

  ~TrackFitterAlgoLTF() {}

  int fit(const TTRoadComb& acomb, TTTrack2& atrack);

 private:
  std::shared_ptr<LinearizedTrackFitter> linearizedTrackFitter_;
  std::shared_ptr<LinearizedTrackFitterEmulator> linearizedTrackFitterEmulator_;
  int verbose_;
  int emu_;
  double normChi2_;

/**
 * Fill the TTTrack2. This could be done with polymorphism on the linearizedTrackFitter and
 * linearizedTrackFitterEmulator. However, the signature of the fit function is different as they take integers
 * or floats for the variables (virtual template functions are not allowed in C++). We prefer to template
 * the common part of the track filling.
 */
  template <class T>
  void fillTrack(T& fitter, TTTrack2& atrack)
  {
    const std::vector<double> &pars = fitter->estimatedPars();
    //const std::vector<double>& principals = fitter->principalComponents();
    const std::vector<double> &principals = fitter->normalizedPrincipalComponents();
    int ndof = fitter->ndof();
    atrack.setTrackParams(0.003 * 3.8 * pars[0], pars[1], pars[2], pars[3], 0, normChi2_ * ndof, ndof, 0, 0);
    std::vector<float> principals_vec;
    for (unsigned ivar = 0; ivar < principals.size(); ++ivar) {
      principals_vec.push_back(principals.at(ivar));
    }
    atrack.setPrincipals(principals_vec);
  }
};

}  // namespace slhcl1tt

#endif
