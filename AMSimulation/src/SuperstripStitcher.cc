#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/SuperstripStitcher.h"
using namespace slhcl1tt;

//#include <algorithm>
#include <iostream>

static const unsigned MAX_NLAYERS = pattern_type().size();


void SuperstripStitcher::init() {
    if (!(6 <= nLayers_ && nLayers_ <= 8))
        throw std::invalid_argument("Incorrect nLayers is given.");

    if (!(nFakeSuperstrips_ <= 3))
        throw std::invalid_argument("Incorrect nFakeSuperstrip is given.");


    unsigned seven_choose_six[7 * 6] = {
        0, 1, 2, 3, 4, 5,
        0, 1, 2, 3, 4, 6,
        0, 1, 2, 3, 5, 6,
        0, 1, 2, 4, 5, 6,
        0, 1, 3, 4, 5, 6,
        0, 2, 3, 4, 5, 6,
        1, 2, 3, 4, 5, 6
    };

    unsigned eight_choose_six[8 * 7 / 2 * 6] = {
        0, 1, 2, 3, 4, 5,
        0, 1, 2, 3, 4, 6,
        0, 1, 2, 3, 4, 7,
        0, 1, 2, 3, 5, 6,
        0, 1, 2, 3, 5, 7,
        0, 1, 2, 3, 6, 7,
        0, 1, 2, 4, 5, 6,
        0, 1, 2, 4, 5, 7,
        0, 1, 2, 4, 6, 7,
        0, 1, 2, 5, 6, 7,
        0, 1, 3, 4, 5, 6,
        0, 1, 3, 4, 5, 7,
        0, 1, 3, 4, 6, 7,
        0, 1, 3, 5, 6, 7,
        0, 1, 4, 5, 6, 7,
        0, 2, 3, 4, 5, 6,
        0, 2, 3, 4, 5, 7,
        0, 2, 3, 4, 6, 7,
        0, 2, 3, 5, 6, 7,
        0, 2, 4, 5, 6, 7,
        0, 3, 4, 5, 6, 7,
        1, 2, 3, 4, 5, 6,
        1, 2, 3, 4, 5, 7,
        1, 2, 3, 4, 6, 7,
        1, 2, 3, 5, 6, 7,
        1, 2, 4, 5, 6, 7,
        1, 3, 4, 5, 6, 7,
        2, 3, 4, 5, 6, 7
    };

    unsigned eight_choose_seven[8 * 7] = {
        0, 1, 2, 3, 4, 5, 6,
        0, 1, 2, 3, 4, 5, 7,
        0, 1, 2, 3, 4, 6, 7,
        0, 1, 2, 3, 5, 6, 7,
        0, 1, 2, 4, 5, 6, 7,
        0, 1, 3, 4, 5, 6, 7,
        0, 2, 3, 4, 5, 6, 7,
        1, 2, 3, 4, 5, 6, 7
    };

    std::map<unsigned, unsigned> layer_map = {
        {5,0}, {6,1}, {7,2}, {8,3}, {9,4}, {10,5},
        {11,5}, {12,4}, {13,3}, {14,2}, {15,1},
        {18,5}, {19,4}, {20,3}, {21,2}, {22,1}
    };

    // Assign
    seven_choose_six_.assign(seven_choose_six, seven_choose_six+(7 * 6));
    eight_choose_six_.assign(eight_choose_six, eight_choose_six+(8 * 7 / 2 * 6));
    eight_choose_seven_.assign(eight_choose_seven, eight_choose_seven+(8 * 7));

    layer_map_.insert(layer_map.begin(), layer_map.end());
}

static std::vector<pattern_type> ret;
static std::vector<pattern_type> ret2;
std::vector<pattern_type> SuperstripStitcher::stitch(const std::vector<addr_type>& superstrips) const {
    ret.clear();

    // nFakeSuperstrips is not enough, return empty vector
    if (nLayers_ >  superstrips.size() + nFakeSuperstrips_)
        return ret;

    // Copy up to 8 superstrips
    ret.resize(1);
    std::copy(superstrips.begin(),
              superstrips.begin() + (MAX_NLAYERS < superstrips.size() ? MAX_NLAYERS : superstrips.size()),
              ret.back().begin() );

    // Case 0: Exact size
    if (nLayers_ == superstrips.size())
        return ret;

    // Case 1: Need fake superstrips
    unsigned i = superstrips.size();
    while (nLayers_ > i) {  // lack of 1 - 3 superstrips
        ret.back().at(i) = fakeSuperstripId_ - (i - superstrips.size());
        i++;
    }
    if (nLayers_ == i)
        return ret;

    if (strategy_ == 0) {
        // Case 2: Do layer merging
        // Prefer endcap to barrel
        ret.back().fill(0);
        id_type lay;
        for (unsigned i=0; i<superstrips.size(); ++i) {
            lay = decodeSuperstripLayer(superstrips.at(i));
            ret.back().at(layer_map_.at(lay)) = superstrips.at(i);
        }
        //std::sort(ret.back().begin(), ret.back().begin() + nLayers_, std::less<addr_type>());
        return ret;

    } else {
        // Case 2: Need to do all combinations
        // The combinations are hardcoded for nLayers = 6, 7, 8
        const std::vector<unsigned>* n_choose_k_ptr = 0;
        if (superstrips.size() == 7 && nLayers_ == 6) {
            n_choose_k_ptr = &seven_choose_six_;
        } else if (superstrips.size() == 8 && nLayers_ == 6) {
            n_choose_k_ptr = &eight_choose_six_;
        } else if (superstrips.size() == 8 && nLayers_ == 7) {
            n_choose_k_ptr = &eight_choose_seven_;
        } else {
            std::cout << "ERROR: Not coded to do '" << superstrips.size() << " choose " << nLayers_ << "' combination. Sorry." << std::endl;
            throw std::logic_error("Unexpected combination of superstrips.");
        }

        ret2.clear();
        ret2.resize(n_choose_k_ptr -> size() / nLayers_);
        unsigned j = 0;
        for (i=0; i<ret2.size(); i++) {
            for (j=0; j<nLayers_; j++) {
                ret2.at(i).at(j) = ret.back().at(n_choose_k_ptr -> at(i * nLayers_ + j));
            }
        }
        return ret2;
    }
}

void SuperstripStitcher::print() {
    std::cout << "nLayers: " << nLayers_ << " nFakeSuperstrips: " << nFakeSuperstrips_ << std::endl;
}
