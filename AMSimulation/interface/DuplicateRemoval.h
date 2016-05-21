#ifndef AMSimulation_DuplicateRemoval_h_
#define AMSimulation_DuplicateRemoval_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"
#include <vector>
#include <algorithm>


namespace slhcl1tt {

class DuplicateRemoval {
  public:
    // Constructor
    DuplicateRemoval() {}

    // Destructor
    ~DuplicateRemoval() {}

    // Return flags categorizing as duplicate (1) or not (0)
    void CheckTracks(std::vector<TTTrack2>& full_am_track_list, int dupRm);
};

}

#endif
