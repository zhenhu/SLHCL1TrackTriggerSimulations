#ifndef AMSimulationDataFormats_TTRoadComb_h_
#define AMSimulationDataFormats_TTRoadComb_h_

#include <vector>
#include <iosfwd>


namespace slhcl1tt {

struct TTRoadComb {
    unsigned roadRef;
    unsigned patternRef;
    unsigned tower;
    unsigned nstubs;
    std::vector<unsigned> stubRefs;

    std::vector<float> stubs_r;
    std::vector<float> stubs_phi;
    std::vector<float> stubs_z;
    std::vector<bool>  stubs_bool;
};


// _____________________________________________________________________________
// Output streams
std::ostream& operator<<(std::ostream& o, const TTRoadComb& road);

}  // namespace slhcl1tt

#endif
