#ifndef AMSimulation_TrackFitterAlgoBase_h_
#define AMSimulation_TrackFitterAlgoBase_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTHit.h"
//#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTTrack.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTTrack2.h"

#include <iosfwd>
#include <string>
#include <vector>


class TrackFitterAlgoBase {
  public:
    TrackFitterAlgoBase() {}
    ~TrackFitterAlgoBase() {}

    //virtual int fit(const std::vector<TTHit>& hits, TTTrack2& track) = 0;

};

#endif

