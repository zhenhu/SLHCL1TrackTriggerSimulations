#ifndef AMSimulation_PatternBankOption_h_
#define AMSimulation_PatternBankOption_h_

#include <string>
#include <vector>
#include <iosfwd>

namespace slhcl1tt {

struct PatternBankOption {
    int         verbose;

    std::string input;
    std::string output;
    std::string bankfile;
    std::string matrixfile;
    std::string roadfile;
    std::string trackfile;

    long long   maxEvents;
    unsigned    nLayers;
    unsigned    nFakers;
    unsigned    nDCBits;

    unsigned    tower;
    std::string superstrip;
    std::string algo;

    float       minPt;
    float       maxPt;
    float       minEta;
    float       maxEta;
    float       minPhi;
    float       maxPhi;
    float       minVz;
    float       maxVz;

    int         minFrequency;
    long int    maxPatterns;
    int         maxMisses;
    int         maxStubs;
    int         maxRoads;

    unsigned    hitbits;

    float       maxChi2;
    int         minNdof;
    int         maxCombs;
    int         maxTracks;

    std::string datadir;
};

std::ostream& operator<<(std::ostream& o, const PatternBankOption& po);

}  // namespace slhcl1tt

#endif

