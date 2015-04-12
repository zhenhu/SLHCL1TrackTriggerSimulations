#ifndef AMSimulation_PatternBankOption_h_
#define AMSimulation_PatternBankOption_h_

#include <string>
#include <vector>
#include <iosfwd>

namespace slhcl1tt {

struct PatternBankOption {
    int verbose;

    std::string input, output, datadir, bankfile, roadfile, trackfile;
    long long maxEvents;
    int minFrequency, maxPatterns, maxMisses, maxStubs, maxRoads, maxCombs, maxTracks;
    bool notrim;

    std::string superstrip;
    float minPt;
    float maxPt;
    float minEta;  // not absolute eta
    float maxEta;  // not absolute eta
    float minPhi;  // from -pi to pi
    float maxPhi;  // from -pi to pi
    unsigned nLayers;
    unsigned nFakers;
    unsigned nDCBits;
    unsigned tower;

    std::string algo;
    float    maxChi2;
    int      minNdof;
};

std::ostream& operator<<(std::ostream& o, const PatternBankOption& po);

}  // namespace slhcl1tt

#endif

