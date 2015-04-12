#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternBankOption.h"

#include <iostream>
#include <iterator>


namespace slhcl1tt {

// Print vector
template<class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(o, " "));
    return o;
}

std::ostream& operator<<(std::ostream& o, const PatternBankOption& po) {
    o << "verbose: "        << po.verbose

      << "  input: "        << po.input
      << "  output: "       << po.output
      << "  bankfile: "     << po.bankfile
      << "  matrixfile: "   << po.matrixfile
      << "  roadfile: "     << po.roadfile
      << "  trackfile: "    << po.trackfile

      << "  maxEvents: "    << po.maxEvents
      << "  nLayers: "      << po.nLayers
      << "  nFakers: "      << po.nFakers
      << "  nDCBits: "      << po.nDCBits

      << "  tower: "        << po.tower
      << "  superstrip: "   << po.superstrip
      << "  algo: "         << po.algo

      << "  minPt: "        << po.minPt
      << "  maxPt: "        << po.maxPt
      << "  minEta: "       << po.minEta
      << "  maxEta: "       << po.maxEta
      << "  minPhi: "       << po.minPhi
      << "  maxPhi: "       << po.maxPhi
      << "  minVz: "        << po.minVz
      << "  maxVz: "        << po.maxVz

      << "  minFrequency: " << po.minFrequency
      << "  maxPatterns: "  << po.maxPatterns
      << "  maxMisses: "    << po.maxMisses
      << "  maxStubs: "     << po.maxStubs
      << "  maxRoads: "     << po.maxRoads

      << "  maxChi2: "      << po.maxChi2
      << "  minNdof: "      << po.minNdof
      << "  maxCombs: "     << po.maxCombs
      << "  maxTracks: "    << po.maxTracks

      << "  datadir: "      << po.datadir
      ;
    return o;
}

}
